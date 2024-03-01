#[derive(Debug, Clone)]
pub struct Region {
    pub seqname: String,
    pub start: i64,
    pub end: i64,
}

pub fn sort_regions_in_place(regions: &mut Vec<Region>) {
    regions.sort_by(|a, b| {
        a.seqname
            .cmp(&b.seqname)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
    });
}

pub fn compress_regions(regions: &Vec<Region>) -> Vec<Region> {
    let mut compressed = vec![];
    let mut current = regions[0].clone();
    for region in regions.iter().skip(1) {
        if region.seqname == current.seqname && region.start <= current.end {
            current.end = region.end;
        } else {
            compressed.push(current.clone());
            current = region.clone();
        }
    }
    compressed.push(current);
    compressed
}