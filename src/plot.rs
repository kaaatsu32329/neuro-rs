use crate::hodgkin_huxley::HodgkinHuxley;
use plotters::prelude::*;

pub trait Plot {
    fn plot(&self, name: &str) -> Result<(), Box<dyn std::error::Error>>;
}

impl Plot for Vec<f32> {
    fn plot(&self, name: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut xs = Vec::new();
        let ys = self;
        for i in 0..self.len() {
            xs.push(i as u32);
        }

        let width = 1980;
        let height = 1080;
        let path = format!("{}{}{}", "./graph/", name, ".png");
        let root = BitMapBackend::new(&path, (width, height)).into_drawing_area();

        root.fill(&WHITE)?;

        let (y_min, y_max) = ys
            .iter()
            .fold((f32::NAN, f32::NAN), |(m, n), v| (v.min(m), v.max(n)));

        let font = ("sans-serif", 20);

        let mut chart;

        if y_min.is_sign_negative() {
            chart = ChartBuilder::on(&root)
                .caption(name, font.into_font())
                .margin(10)
                .x_label_area_size(16)
                .y_label_area_size(42)
                .build_cartesian_2d(*xs.first().unwrap()..*xs.last().unwrap(), y_min..y_max)?;
        } else {
            chart = ChartBuilder::on(&root)
                .caption(name, font.into_font())
                .margin(10)
                .x_label_area_size(16)
                .y_label_area_size(42)
                .build_cartesian_2d(*xs.first().unwrap()..*xs.last().unwrap(), 0f32..y_max)?;
            // ToDo: Update
        }

        chart.configure_mesh().draw()?;

        let line_series = LineSeries::new(xs.iter().zip(ys.iter()).map(|(x, y)| (*x, *y)), &RED);
        chart.draw_series(line_series)?;

        Ok(())
    }
}

impl Plot for HodgkinHuxley {
    fn plot(&self, _name: &str) -> Result<(), Box<dyn std::error::Error>> {
        todo!()
    }
}
