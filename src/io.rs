use crate::regions::{sort_regions_in_place, Region};
use csv::Reader;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use anyhow::Error;

pub struct GtfFile {
    pub path: PathBuf,
}

impl GtfFile {
    pub fn new(file_path: impl Into<PathBuf>) -> Self {
        GtfFile {
            path: file_path.into(),
        }
    }

    pub fn reader(&self) -> Result<Reader<BufReader<MultiGzDecoder<File>>>, Error> {
        let file = File::open(&self.path)?;
        let decoder = MultiGzDecoder::new(file);
        let reader = BufReader::new(decoder);
        let csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);
        Ok(csv_reader)
    }

    // Selects regions marked as "exon", transforms their coordinates into
    // 0-based, half-open intervals, sorts them by chromosome and position,
    // and returns them as a vector of Region structs.
    pub fn exon_regions(&self) -> Result<Vec<Region>, Error> {
        let mut regions = vec![];
        let mut reader = self.reader()?;
        for result in reader.records() {
            let record = result?;
            if &record[2] == "exon" {
                regions.push(Region {
                    seqname: record[0].to_string(),
                    start: (record[3].parse::<i64>()?) - 1i64,
                    end: record[4].parse()?,
                });
            }
        }
        sort_regions_in_place(&mut regions);
        Ok(regions)
    }
}
