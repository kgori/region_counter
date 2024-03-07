use clap::error::ErrorKind;
use clap::{CommandFactory, Parser};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug, Clone)]
#[command(version, about, long_about = None)]
pub struct ProgramOptions {
    #[arg(short, long)]
    pub bamfile: PathBuf,

    #[arg(short = 'g', long)]
    pub gtf: PathBuf,

    #[arg(short = 'q', long, default_value = "35")]
    pub minmapqual: u8,

    #[arg(short = 'f', long, default_value = "3")]
    pub required_flag: u16,

    #[arg(short = 'F', long, default_value = "2816")]
    pub filtered_flag: u16,
}

fn validate_file(file: &Path) {
    if !file.exists() {
        let mut cmd = ProgramOptions::command();
        cmd.error(
            ErrorKind::ValueValidation,
            format!("file `{}` not found", file.display()),
        )
        .exit();
    }
}

pub fn parse_cli() -> ProgramOptions {
    let args = ProgramOptions::parse();
    validate_file(&args.bamfile);
    validate_file(&args.gtf);
    args
}
