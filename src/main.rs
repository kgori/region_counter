#![allow(unused_imports)]
#![allow(dead_code)]
use cigar::check_cigar_overlap;
use cli::ProgramOptions;
use pdatastructs::hyperloglog::HyperLogLog;
use regions::{compress_regions, convert_regions_vec_to_hashmap, Region};
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use rayon::prelude::*;
use anyhow::Error;

const ADDRESSBITS: usize = 18;

mod cigar;
mod cli;
mod io;
mod regions;

fn get_read_name(read: &Record) -> String {
    let mut read_name = std::str::from_utf8(read.qname()).unwrap().to_owned();
    if read.is_first_in_template() {
        read_name.push_str("/1");
    } else if read.is_last_in_template() {
        read_name.push_str("/2");
    }
    read_name
}

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
    let chroms = chroms.iter().map(|x| String::from_utf8(x.to_vec())).collect::<Result<Vec<_>, _>>();
    match chroms {
        Ok(chroms) => Ok(chroms),
        Err(e) => Err(Error::msg(format!("Error reading chromosome names: {}", e))),
    }
}

fn count_reads(
    chrom: &str,
    regions: &Vec<Region>,
    args: &ProgramOptions,
) -> Result<(usize, usize), Error> {
    // Only count unique reads
    let mut all_reads_set = HashSet::new();
    let mut exon_reads_set = HashSet::new();

    let mut bam = IndexedReader::from_path(&args.bamfile)?;
    let mut read = Record::new();
    bam.fetch(chrom)?;

    let mut current_region_index = 0;
    let max_index = regions.len();

    while let Some(result) = bam.read(&mut read) {
        match result {
            Ok(_) => {
                match check_read(&read, args) {
                    ReadCheckOutcome::Reject => continue,
                    ReadCheckOutcome::Accept => {},
                }

                let read_name = get_read_name(&read);
                all_reads_set.insert(read_name.clone());

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
                            exon_reads_set.insert(read_name);
                            break;
                        }
                    }
                }
            }
            Err(e) => println!("Error reading read: {}", e),
        }
    }
    Ok((all_reads_set.len(), exon_reads_set.len()))
}

fn count_mapped_reads(
    args: &ProgramOptions,
    regions: &HashMap<String, Vec<Region>>,
) -> Result<(usize, usize), Error> {
    let mut all_reads = 0;
    let mut exon_reads = 0;
    
    let mut chroms: Vec<_> = regions.keys().collect();
    chroms.sort();

    let results: Vec<Result<(usize, usize), Error>> = chroms.par_iter().map(|chrom| {
        println!("Counting reads on chromosome {}", chrom);
        let regions = regions.get(*chrom).unwrap();
        count_reads(chrom, regions, args)
    }).collect();

    for result in results {
        let (all, exon) = result?;
        all_reads += all;
        exon_reads += exon;
    }

    Ok((all_reads, exon_reads))
}

fn count_unmapped_reads(args: &ProgramOptions) -> Result<usize, Error> {
    let mut bam = IndexedReader::from_path(&args.bamfile).unwrap();
    bam.fetch("*")?;
    let mut read = Record::new();
    let mut unmapped = HashSet::new();
    while let Some(result) = bam.read(&mut read) {
        match result {
            Ok(_) => {
                if read.is_unmapped() && !read.is_quality_check_failed() && !read.is_duplicate() {
                    let read_name = get_read_name(&read);
                    unmapped.insert(read_name);
                }
            }
            Err(e) => println!("Error reading read: {}", e),
        }
    }
    Ok(unmapped.len())
}

fn main() -> Result<(), Error> {
    let args = cli::parse_cli();
    let gtf = io::GtfFile::new(&args.gtf);
    println!("Reading GTF file: {}", args.gtf.display());
    let regions = gtf.exon_regions()?;
    let regions = compress_regions(&regions);
    let mut regions_map = convert_regions_vec_to_hashmap(regions);
    let n_regions = regions_map
        .iter()
        .map(|(_, v)| {
            v.len()
        })
        .sum::<usize>();
    let chroms = get_chrom_names(&args.bamfile)?;
    for chrom in chroms {
        if !regions_map.contains_key(&chrom) {
            regions_map.insert(chrom, Vec::new());
        }
    }
    println!("Counting {} exon regions on {} chromosomes", n_regions, regions_map.len());
    let (all_reads, exon_reads) = count_mapped_reads(&args, &regions_map)?;
    let unmapped_reads = count_unmapped_reads(&args)?;
    println!("{} total mapped reads", all_reads);
    println!("{} exon mapped reads", exon_reads);
    println!("{} unmapped reads", unmapped_reads);
    Ok(())
}
