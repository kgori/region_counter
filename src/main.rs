#![allow(unused_imports)]
use cli::ProgramOptions;
use regions::compress_regions;
use pdatastructs::hyperloglog::HyperLogLog;
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use std::collections::HashSet;

const ADDRESSBITS: usize = 18;

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

fn count_whole_genome(args: &ProgramOptions) {
    let mut hll = HyperLogLog::new(ADDRESSBITS);
    let mut set = HashSet::new();

    let mut bam = Reader::from_path(&args.bamfile).unwrap();
    let mut read = Record::new();
    let mut nreads: u32 = 0;
    while let Some(result) = bam.read(&mut read) {
        match result {
            Ok(_) => {
                let read_name = get_read_name(&read);

                match check_read(&read, args) {
                    ReadCheckOutcome::Reject => continue,
                    ReadCheckOutcome::Accept => {
                        nreads += 1;
                        hll.add(&read_name);
                        set.insert(read_name.clone());
                    }
                }
            }
            Err(e) => println!("Error reading read: {}", e),
        }
    }
    println!("{} reads", nreads);
    println!("{} unique reads", set.len());
    println!("{} unique reads (HLL)", hll.count());
}

fn count_exons(args: &ProgramOptions, regions: &Vec<regions::Region>) {
    let mut hll = HyperLogLog::new(ADDRESSBITS);
    let mut set = HashSet::new();

    let mut bam = IndexedReader::from_path(&args.bamfile).unwrap();
    let mut read = Record::new();
    let mut nreads: u32 = 0;

    for region in regions {
        println!("Processing region: {:?}", region);
        bam.fetch((&region.seqname, region.start, region.end))
            .unwrap();
        while let Some(result) = bam.read(&mut read) {
            match result {
                Ok(_) => {
                    let read_name = get_read_name(&read);

                    match check_read(&read, args) {
                        ReadCheckOutcome::Reject => continue,
                        ReadCheckOutcome::Accept => {
                            nreads += 1;
                            hll.add(&read_name);
                            set.insert(read_name.clone());
                        }
                    }
                }
                Err(e) => println!("Error reading read: {}", e),
            }
        }
    }

    println!("{} reads", nreads);
    println!("{} unique reads", set.len());
    println!("{} unique reads (HLL)", hll.count());
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = cli::parse_cli();
    let gtf = io::GtfFile::new(&args.gtf);
    let regions = gtf.exon_regions()?;
    let regions = compress_regions(&regions);
    count_exons(&args, &regions);
    // count_whole_genome(&args);
    Ok(())
}
