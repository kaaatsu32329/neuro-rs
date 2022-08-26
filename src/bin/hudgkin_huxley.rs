use neuro_rs::neuro::HodgkinHuxley;
use neuro_rs::plot::Plot;

fn main() {
    let mut hh_model = HodgkinHuxley::new();
    hh_model.calc();
    // println!("{:?}", hh_model);

    hh_model.voltage.plot("voltage").unwrap();
    hh_model.i_na.plot("Na+").unwrap();
    // ToDo: Update
}
