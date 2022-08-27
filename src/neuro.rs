// Reference: https://github.com/kaaatsu32329/hodgkin-huxley-sim

const CONDUCTANCE_NA: f32 = 120.0;
const CONDUCTANCE_K: f32 = 36.0;
const CONDUCTANCE_LEAK: f32 = 0.3;

const EQUILIBRIUM_VOLTAGE_NA: f32 = 50.0;
const EQUILIBRIUM_VOLTAGE_K: f32 = -77.0;
const EQUILIBRIUM_VOLTAGE_LEAK: f32 = -54.387;

const CAPACITANCE: f32 = 1.0;

const STEP: usize = 3;

#[derive(Debug, Clone, Default)]
pub struct HodgkinHuxley {
    time: usize, // Experimental time
    pub dt: f32, // time step
    current_step: usize,
    pub voltage: Vec<f32>,
    pub i_na: Vec<f32>,
    pub i_k: Vec<f32>,
    pub i_leak: Vec<f32>,
    pub i_injection: Vec<f32>,
    pub m: Vec<f32>,
    pub h: Vec<f32>,
    pub n: Vec<f32>,
}

impl HodgkinHuxley {
    pub fn new() -> Self {
        Self {
            time: 200000,
            dt: 0.001,
            current_step: 0,
            voltage: vec![-65.],
            i_na: vec![0.],
            i_k: vec![0.],
            i_leak: vec![0.],
            i_injection: vec![0.],
            m: vec![0.05],
            h: vec![0.6],
            n: vec![0.32],
        }
    }

    pub fn calc(&mut self) {
        for _ in 0..self.time {
            self.voltage
                .push(self.voltage[self.current_step] + self.diff_membrane_potential() * self.dt);
            self.m.push(
                self.m[self.current_step]
                    + self.diff_gate(self.alpha_m(), self.beta_m(), self.m[self.current_step])
                        * self.dt,
            );
            self.h.push(
                self.h[self.current_step]
                    + self.diff_gate(self.alpha_h(), self.beta_h(), self.h[self.current_step])
                        * self.dt,
            );
            self.n.push(
                self.n[self.current_step]
                    + self.diff_gate(self.alpha_n(), self.beta_n(), self.n[self.current_step])
                        * self.dt,
            );
            self.intensity_injection();
            self.intensity_na();
            self.intensity_k();
            self.intensity_leak();
            self.current_step += 1;

            #[cfg(feature="print")]
            {
                println!("voltage: {}", self.voltage[self.current_step]);
                println!(
                    "i: {}, {}, {}",
                    self.i_na[self.current_step],
                    self.i_k[self.current_step],
                    self.i_leak[self.current_step]
                );
                println!(
                    "gate: {}, {}, {}\n",
                    self.m[self.current_step], self.h[self.current_step], self.n[self.current_step]
                );
            }
        }
    }

    fn intensity_na(&mut self) {
        self.i_na.push(
            CONDUCTANCE_NA
                * self.m[self.current_step].powi(3)
                * self.h[self.current_step]
                * self.voltage[self.current_step]
                - EQUILIBRIUM_VOLTAGE_NA,
        );
    }

    fn intensity_k(&mut self) {
        self.i_k.push(
            CONDUCTANCE_K
                * self.n[self.current_step].powi(4)
                * (self.voltage[self.current_step] - EQUILIBRIUM_VOLTAGE_K),
        )
    }

    fn intensity_leak(&mut self) {
        self.i_leak
            .push(CONDUCTANCE_LEAK * (self.voltage[self.current_step] - EQUILIBRIUM_VOLTAGE_LEAK))
    }

    fn intensity_injection(&mut self) {
        let y = (self.dt * (self.current_step * (STEP + 1)) as f32 / 5.).sin();
        if y.is_sign_positive() {
            self.i_injection.push(30.);
        } else {
            self.i_injection.push(0.);
        }
    }

    fn alpha_m(&self) -> f32 {
        0.1 * (-40. - self.voltage[self.current_step])
            / (((-40. - self.voltage[self.current_step]) / 10.).exp() - 1.)
    }

    fn alpha_h(&self) -> f32 {
        0.07 * ((-self.voltage[self.current_step] - 65.) / 20.).exp()
    }

    fn alpha_n(&self) -> f32 {
        0.01 * (-55. - self.voltage[self.current_step])
            / (((-55. - self.voltage[self.current_step]) / 10.).exp() - 1.)
    }

    fn beta_m(&self) -> f32 {
        4. * ((-self.voltage[self.current_step] - 65.) / 18.).exp()
    }

    fn beta_h(&self) -> f32 {
        1. / (((-35. - self.voltage[self.current_step]) / 10.).exp() + 1.)
    }

    fn beta_n(&self) -> f32 {
        0.125 * ((-65. - self.voltage[self.current_step]) / 80.).exp()
    }

    fn diff_membrane_potential(&self) -> f32 {
        (self.i_injection[self.current_step]
            - self.i_na[self.current_step]
            - self.i_k[self.current_step]
            - self.i_leak[self.current_step])
            / CAPACITANCE
    }

    fn diff_gate(&self, alpha: f32, beta: f32, gate: f32) -> f32 {
        alpha * (1.0 - gate) - beta * gate
    }
}
