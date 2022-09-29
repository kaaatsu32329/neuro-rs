use neuro_rs::hodgkin_huxley::HodgkinHuxley;
use neuro_rs::plot::Plot;

fn main() {
    let mut hh_model = HodgkinHuxley::new();
    hh_model.calc();

    hh_model.voltage.plot("voltage").unwrap();
    hh_model.i_na.plot("Na+").unwrap();
    // ToDo: Update
}
