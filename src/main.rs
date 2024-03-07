use anyhow::Error;
use cigar::check_cigar_overlap;
use cli::ProgramOptions;
use rayon::prelude::*;
use regions::{compress_regions, convert_regions_vec_to_hashmap, Region};
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::collections::HashMap;

mod cigar;
mod cli;
mod io;
mod regions;

const FLAGS_ALWAYS_FILTERED: u16 = 2816;
const FLAG_PROPER_PAIR: u16 = 2;
const FLAG_UNMAPPED: u16 = 4;
const FLAGS_MAPPING_RELATED: u16 = 63;

enum ReadCheckOutcome {
    Accept,
    Reject,
}

fn check_read(read: &Record, args: &ProgramOptions) -> ReadCheckOutcome {
    if read.mapq() < args.minmapqual {
        return ReadCheckOutcome::Reject;
    }

    if read.flags() & args.required_flag != args.required_flag
        || read.flags() & args.filtered_flag != 0
    {
        return ReadCheckOutcome::Reject;
    }

    ReadCheckOutcome::Accept
}

fn get_chrom_names(bamfile: &std::path::Path) -> Result<Vec<String>, Error> {
    let bam = Reader::from_path(bamfile)?;
    let header = bam.header();
    let chroms = header.target_names();
    let chroms = chroms
        .iter()
        .map(|x| String::from_utf8(x.to_vec()))
        .collect::<Result<Vec<_>, _>>();
    match chroms {
        Ok(chroms) => Ok(chroms),
        Err(e) => Err(Error::msg(format!("Error reading chromosome names: {}", e))),
    }
}

struct CountResult {
    accepted: usize,
    rejected: usize,
}

fn count_reads(
    chrom: &str,
    regions: &Vec<Region>,
    args: &ProgramOptions,
) -> Result<(CountResult, CountResult), Error> {
    // Only count unique reads
    let mut all_reads_accepted = 0;
    let mut exon_reads_accepted = 0;
    let mut all_reads_rejected = 0;
    let mut exon_reads_rejected = 0;

    let mut bam = IndexedReader::from_path(&args.bamfile)?;
    let mut read = Record::new();
    bam.fetch(chrom)?;

    let mut current_region_index = 0;
    let max_index = regions.len();

    while let Some(result) = bam.read(&mut read) {
        match result {
            Ok(_) => {
                if read.flags() & FLAGS_ALWAYS_FILTERED != 0 {
                    continue;
                }

                let read_check_outcome = check_read(&read, args);

                match read_check_outcome {
                    ReadCheckOutcome::Accept => {
                        all_reads_accepted += 1;
                    }
                    ReadCheckOutcome::Reject => {
                        all_reads_rejected += 1;
                    }
                }

                // Check if the read is past the end of the current region
                // If it is, advance to the next region as long as there are regions left
                loop {
                    if current_region_index >= max_index {
                        break;
                    }
                    let region = &regions[current_region_index];
                    if read.pos() >= region.end {
                        current_region_index += 1;
                    } else {
                        break;
                    }
                }

                // If there still is a current region, check if the read overlaps it
                if current_region_index < max_index {
                    let end_pos = cigar::cigar_end_pos(&read);
                    for index in current_region_index..max_index {
                        let region = &regions[index];
                        if region.start > end_pos {
                            break;
                        }
                        if check_cigar_overlap(&read, region.start, region.end) {
                            match read_check_outcome {
                                ReadCheckOutcome::Accept => {
                                    exon_reads_accepted += 1;
                                }
                                ReadCheckOutcome::Reject => {
                                    exon_reads_rejected += 1;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            Err(e) => println!("Error reading read: {}", e),
        }
    }
    let mapped_result = CountResult {
        accepted: all_reads_accepted,
        rejected: all_reads_rejected,
    };
    let exon_result = CountResult {
        accepted: exon_reads_accepted,
        rejected: exon_reads_rejected,
    };
    Ok((mapped_result, exon_result))
}

fn count_mapped_reads(
    args: &ProgramOptions,
    regions: &HashMap<String, Vec<Region>>,
) -> Result<(CountResult, CountResult), Error> {
    let mut all_reads = CountResult {
        accepted: 0,
        rejected: 0,
    };
    let mut exon_reads = CountResult {
        accepted: 0,
        rejected: 0,
    };

    let mut chroms: Vec<_> = regions.keys().collect();
    chroms.sort();

    let results: Vec<Result<(CountResult, CountResult), Error>> = chroms
        .par_iter()
        .map(|chrom| {
            eprintln!("Counting reads on chromosome {}", chrom);
            let regions = regions.get(*chrom).unwrap();
            count_reads(chrom, regions, args)
        })
        .collect();

    for result in results {
        let (all, exon) = result?;
        all_reads.accepted += all.accepted;
        all_reads.rejected += all.rejected;
        exon_reads.accepted += exon.accepted;
        exon_reads.rejected += exon.rejected;
    }

    Ok((all_reads, exon_reads))
}

fn count_unmapped_reads(args: &ProgramOptions) -> Result<CountResult, Error> {
    let mut args: ProgramOptions = args.clone();
    args.minmapqual = 0;
    args.required_flag ^= args.required_flag & FLAG_PROPER_PAIR; // Turn off mapping requirement
    args.required_flag ^= FLAG_UNMAPPED; // Turn on unmapped requirement
    args.filtered_flag ^= args.filtered_flag & FLAGS_MAPPING_RELATED; // Turn off mapping related flags
    let mut bam = IndexedReader::from_path(&args.bamfile).unwrap();
    bam.fetch("*")?;
    let mut read = Record::new();
    let mut unmapped_accepted = 0;
    let mut unmapped_rejected = 0;
    while let Some(result) = bam.read(&mut read) {
        match result {
            Ok(_) => {
                if read.flags() & FLAGS_ALWAYS_FILTERED != 0 {
                    continue;
                }

                let read_check_outcome = check_read(&read, &args);
                match read_check_outcome {
                    ReadCheckOutcome::Accept => {
                        unmapped_accepted += 1;
                    }
                    ReadCheckOutcome::Reject => {
                        unmapped_rejected += 1;
                    }
                }
            }
            Err(e) => println!("Error reading read: {}", e),
        }
    }
    Ok(CountResult {
        accepted: unmapped_accepted,
        rejected: unmapped_rejected,
    })
}

fn main() -> Result<(), Error> {
    let args = cli::parse_cli();
    let gtf = io::GtfFile::new(&args.gtf);
    eprintln!("Reading GTF file: {}", args.gtf.display());
    let regions = gtf.exon_regions()?;
    let regions = compress_regions(&regions);
    let mut regions_map = convert_regions_vec_to_hashmap(regions);
    let n_regions = regions_map.iter().map(|(_, v)| v.len()).sum::<usize>();
    let chroms = get_chrom_names(&args.bamfile)?;
    for chrom in chroms {
        if !regions_map.contains_key(&chrom) {
            regions_map.insert(chrom, Vec::new());
        }
    }
    eprintln!(
        "Counting {} exon regions on {} chromosomes",
        n_regions,
        regions_map.len()
    );
    let (all_reads, exon_reads) = count_mapped_reads(&args, &regions_map)?;
    let unmapped_reads = count_unmapped_reads(&args)?;
    println!("## Min mapping quality: {}", args.minmapqual);
    println!("## Required flag: {}", args.required_flag);
    println!("## Filtered flag: {}", args.filtered_flag);
    println!("## GTF file: {}", args.gtf.display());
    println!("## BAM file: {}", args.bamfile.display());
    println!("Category\tAccepted\tRejected\tTotal");
    println!(
        "Exon\t{}\t{}\t{}",
        exon_reads.accepted,
        exon_reads.rejected,
        exon_reads.accepted + exon_reads.rejected
    );
    println!(
        "Mapped\t{}\t{}\t{}",
        all_reads.accepted,
        all_reads.rejected,
        all_reads.accepted + all_reads.rejected
    );
    println!(
        "Unmapped\t{}\t{}\t{}",
        unmapped_reads.accepted,
        unmapped_reads.rejected,
        unmapped_reads.accepted + unmapped_reads.rejected
    );
    println!(
        "Total\t{}\t{}\t{}",
        all_reads.accepted + unmapped_reads.accepted,
        all_reads.rejected + unmapped_reads.rejected,
        all_reads.accepted + all_reads.rejected + unmapped_reads.accepted + unmapped_reads.rejected
    );
    Ok(())
}
