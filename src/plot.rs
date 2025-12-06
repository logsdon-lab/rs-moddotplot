use std::{f32, ops::Range};

use itertools::Itertools;
use plotters::{
    coord::{types::RangedCoordf32, Shift},
    prelude::*,
};

use crate::Row;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Label {
    pub color: String,
    pub label: String,
}

fn init_chart<'a>(
    ax: &DrawingArea<BitMapBackend<'a>, Shift>,
    title: Option<&str>,
    xlim: Range<f32>,
    ylim: Range<f32>,
    hide_x: bool,
    hide_y: bool,
) -> ChartContext<'a, BitMapBackend<'a>, Cartesian2d<RangedCoordf32, RangedCoordf32>> {
    let mut binding = ChartBuilder::on(ax);
    let chart_avg = binding
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 60);

    if let Some(title) = title {
        chart_avg.caption(title, ("sans-serif", 40));
    }
    let mut chart_avg = chart_avg.build_cartesian_2d(xlim, ylim).unwrap();

    if hide_x && hide_y {
    } else if hide_x {
        chart_avg
            .configure_mesh()
            .disable_mesh()
            .disable_x_axis()
            .draw()
            .unwrap();
    } else if hide_y {
        chart_avg
            .configure_mesh()
            .disable_mesh()
            .disable_y_axis()
            .draw()
            .unwrap();
    } else {
        chart_avg.configure_mesh().disable_mesh().draw().unwrap();
    }

    chart_avg
}

const IDENT_BRKPTS: [(f32, RGBColor); 12] = [
    (90.0, RGBColor(75, 57, 145)),    // "#4b3991"
    (97.5, RGBColor(41, 116, 175)),   // "#2974af"
    (97.75, RGBColor(74, 157, 168)),  // "#4a9da8"
    (98.0, RGBColor(87, 184, 148)),   // "#57b894"
    (98.25, RGBColor(157, 216, 147)), // "#9dd893"
    (98.5, RGBColor(225, 246, 134)),  // "#e1f686"
    (98.75, RGBColor(255, 255, 178)), // "#ffffb2"
    (99.0, RGBColor(253, 218, 121)),  // "#fdda79"
    (99.25, RGBColor(251, 158, 79)),  // "#fb9e4f"
    (99.5, RGBColor(238, 86, 52)),    // "#ee5634"
    (99.75, RGBColor(201, 39, 62)),   // "#c9273e"
    (100.0, RGBColor(138, 0, 51)),    // "#8a0033"
];

struct IdentityPolygon {
    coords: [(f32, f32); 4],
    color: RGBColor,
    ident_range: (f32, f32),
}

impl IdentityPolygon {
    fn new(coords: [(f32, f32); 4], identity: f32) -> Self {
        let (Ok(idx) | Err(idx)) = IDENT_BRKPTS.binary_search_by(|p| p.0.total_cmp(&identity));
        let (brkpt_1, color) = IDENT_BRKPTS[idx];
        // Default to 0.0 identity
        let brkpt_0 = idx
            .checked_sub(1)
            .map(|idx| IDENT_BRKPTS[idx].0)
            .unwrap_or_default();
        let ident_range = (brkpt_0, brkpt_1);
        IdentityPolygon {
            coords,
            color,
            ident_range,
        }
    }
}
const HT_TRI: u32 = 600;
const TRI_SIDE: f32 = f32::consts::SQRT_2 / 2.0;
const POLY_X: [f32; 4] = [TRI_SIDE, 0.0, -TRI_SIDE, 0.0];
const POLY_Y: [f32; 4] = [0.0, TRI_SIDE, 0.0, -TRI_SIDE];

/// Plot self-identity triangle for one chromosome.
pub fn plot_self_ident_tri<'a, 'b>(
    bedpe: &'a [Row],
    // annot_itvs: impl Iterator<Item = (usize, usize, &'b str)>,
    title: &str,
    outfile: &str,
    invert: bool,
) {
    let window = bedpe
        .first()
        .map(|row| row.query_end - row.query_start)
        .expect("No rows.");
    let mut max_x = f32::MIN;
    let mut max_query_start = usize::MIN;
    let mut coords = Vec::with_capacity(bedpe.len());
    for row in bedpe {
        let first_pos = (row.query_start / window) as f32;
        let second_pos = (row.reference_start / window) as f32;
        let x = first_pos + second_pos;
        let y = -first_pos + second_pos;

        max_x = max_x.max(x);
        max_query_start = max_query_start.max(row.query_start);
        coords.push((x, y, row.perc_id_by_events));
    }
    let scale = max_query_start as f32 / max_x;
    let window = window as f32 / scale;
    let invert = if invert { -1.0 } else { 1.0 };

    let mut poly_coords: Vec<IdentityPolygon> = Vec::with_capacity(bedpe.len() * 4);
    let (mut min_x, mut max_x) = (f32::MAX, f32::MIN);
    let (mut min_y, mut max_y) = (f32::MAX, f32::MIN);
    for (x, y, ident) in coords {
        // I could use zip and collect_array but not sure if faster to just manually unroll
        let x_0 = ((POLY_X[0] * window) + x) * scale;
        let x_1 = ((POLY_X[1] * window) + x) * scale;
        let x_2 = ((POLY_X[2] * window) + x) * scale;
        let x_3 = ((POLY_X[3] * window) + x) * scale;

        let y_0 = ((POLY_Y[0] * window) + y) * window * invert;
        let y_1 = ((POLY_Y[1] * window) + y) * window * invert;
        let y_2 = ((POLY_Y[2] * window) + y) * window * invert;
        let y_3 = ((POLY_Y[3] * window) + y) * window * invert;

        let coords = [(x_0, y_0), (x_1, y_1), (x_2, y_2), (x_3, y_3)];
        for (x, y) in coords {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
        poly_coords.push(IdentityPolygon::new(coords, ident))
    }

    let root = BitMapBackend::new(outfile, (1200, HT_TRI + 100)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let (ax_tri, _ax_annot) = root.split_vertically(HT_TRI);

    let mut chart_tri = init_chart(&ax_tri, Some(title), min_x..max_x, min_y..max_y, true, true);

    // let mut chart_avg = init_chart(&axes[1], None, xlim.clone(), 0.0f32..100.0f32, false, false);

    // Draw by identity group and plot by increasing identity
    for ((group, color), polygons) in &poly_coords
        .into_iter()
        .sorted_unstable_by(|a, b| a.ident_range.0.total_cmp(&b.ident_range.1))
        .chunk_by(|a| (a.ident_range, a.color))
    {
        let label = format!("{}-{}", group.0, group.1);
        // Draw and add legend.
        chart_tri
            .draw_series(
                polygons
                    .into_iter()
                    .map(|polygon| Polygon::new(polygon.coords, color)),
            )
            .unwrap()
            .label(label)
            .legend(move |(x, y)| Rectangle::new([(x - 5, y - 5), (x + 5, y + 5)], color));
    }

    chart_tri
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .margin(10)
        .legend_area_size(10)
        .draw()
        .unwrap();

    // chart_annot
    //     .configure_series_labels()
    //     .position(SeriesLabelPosition::UpperRight)
    //     .margin(10)
    //     .legend_area_size(10)
    //     .draw()
    //     .unwrap();

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().unwrap();
}

#[cfg(test)]
mod test {
    use std::fs;

    use crate::{compute_self_identity, plot::plot_self_ident_tri};

    #[test]
    fn test_draw_tri() {
        let rows = compute_self_identity(
            "data/HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509.fa",
            None,
        );
        let outfile = "HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509.png";
        plot_self_ident_tri(
            &rows,
            "HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509",
            outfile,
            false,
        );
        fs::remove_file(outfile).unwrap();
    }
}
