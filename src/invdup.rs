use std::collections::VecDeque;

use indexmap::{IndexMap, IndexSet};

use crate::Row;

#[derive(Debug)]
pub enum Event {
    Duplication,
    Inversion,
}

impl Event {
    /// Expected change in coordinates in ANI matrix
    #[inline]
    fn dxy(&self) -> (isize, isize) {
        match self {
            Event::Duplication => (-1, 1),
            Event::Inversion => (1, -1),
        }
    }

    fn interval_summary(
        &self,
        event_ident_itvs: &[(usize, usize, f32)],
    ) -> ((usize, usize), (usize, usize), f32) {
        match self {
            Event::Duplication => todo!("{event_ident_itvs:?}"),
            Event::Inversion => todo!("{event_ident_itvs:?}"),
        }
    }
}

#[derive(Debug)]
pub struct EventInterval {
    pub chrom: String,
    pub itv_1: (usize, usize),
    pub itv_2: (usize, usize),
    pub event: Event,
    pub avg_ident: f32,
}

impl EventInterval {
    pub fn header() -> &'static str {
        "#chrom\tstart\tend\tevent"
    }
}

pub fn find_simple_inversion_duplications(rows: &[Row]) -> Vec<EventInterval> {
    let Some((chrom, window)) = rows
        .first()
        .map(|row| (&row.reference_name, row.query_end - row.query_start))
    else {
        return vec![];
    };

    let mut events = vec![];
    let mut aln_mtx: IndexMap<usize, IndexMap<usize, f32>> = IndexMap::default();
    for line in rows {
        let x = line.query_start / window;
        let y = line.reference_start / window;
        let ident = line.perc_id_by_events;
        // Convert position to indices.
        aln_mtx
            .entry(x)
            .and_modify(|rec| {
                rec.insert(y, ident);
            })
            .or_insert_with(|| IndexMap::from_iter([(y, ident)]));
    }
    // 7               +
    // 6             +
    // 5           +   *
    // 4         +   *
    // 3       +   *   *
    // 2     +
    // 1   +   *
    // 0 +       *
    //   0 1 2 3 4 5 6 7
    for x in aln_mtx.keys() {
        // Skip self-diagonal, one window above it, and non-identical regions
        // (x, y, dx, dy)
        let mut positions: VecDeque<(isize, isize, Event)> = VecDeque::from_iter(
            aln_mtx[x]
                .iter()
                .filter(|(y, ident)| *y != x && **y != *x + 1 && **ident != 0.0)
                .flat_map(|(y, _ident)| {
                    [
                        // Move up or sideways
                        (*x as isize, *y as isize, Event::Duplication),
                        (*x as isize, *y as isize, Event::Inversion),
                    ]
                }),
        );
        // Store identity for given event
        let mut event_ident_itvs: Vec<(usize, usize, f32)> = vec![];
        while let Some(position) = positions.pop_front() {
            let (x, y, event) = position;
            let (dx, dy) = event.dxy();

            let adj_x = (x + dx).clamp(0, isize::MAX) as usize;
            let adj_y = (y + dy).clamp(0, isize::MAX) as usize;
            eprintln!("({x}, {y}) -> ({adj_x}, {adj_y})");

            match aln_mtx.get(&adj_x).and_then(|ys| ys.get(&adj_y)) {
                Some(ident) => {
                    // Push front to ensure finish
                    positions.push_front((adj_x as isize, adj_y as isize, event));
                    event_ident_itvs.push((x as usize, y as usize, *ident));
                }
                None => {
                    if event_ident_itvs.is_empty() {
                        continue;
                    }
                    eprintln!("{event:?}: {event_ident_itvs:?}");
                    // let (itv_1, itv_2, avg_ident) = event.interval_summary(&event_ident_itvs);
                    event_ident_itvs.clear();
                    // events.push(EventInterval {
                    //     chrom: chrom.to_owned(),
                    //     itv_1,
                    //     itv_2,
                    //     event,
                    //     avg_ident,
                    // });
                }
            }
        }
    }

    events
}
