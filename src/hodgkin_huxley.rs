// Reference: https://github.com/kaaatsu32329/hodgkin-huxley-sim

const CONDUCTANCE_NA: f32 = 120.0;
const CONDUCTANCE_K: f32 = 36.0;
const CONDUCTANCE_LEAK: f32 = 0.3;

const EQUILIBRIUM_VOLTAGE_NA: f32 = 50.0;
const EQUILIBRIUM_VOLTAGE_K: f32 = -77.0;
const EQUILIBRIUM_VOLTAGE_LEAK: f32 = -54.387;

const CAPACITANCE: f32 = 1.0;

const STEP: usize = 3;

const VOLTAGE_INIT: f32 = -65.;
const GATE_M_INIT: f32 = 0.05;
const GATE_H_INIT: f32 = 0.6;
const GATE_N_INIT: f32 = 0.32;

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
            voltage: vec![VOLTAGE_INIT],
            i_na: vec![0.],
            i_k: vec![0.],
            i_leak: vec![0.],
            i_injection: vec![0.],
            m: vec![GATE_M_INIT],
            h: vec![GATE_H_INIT],
            n: vec![GATE_N_INIT],
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

            #[cfg(feature = "print")]
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
        0.07 * ((-self.voltage[self.current_step] + VOLTAGE_INIT) / 20.).exp()
    }

    fn alpha_n(&self) -> f32 {
        0.01 * (-55. - self.voltage[self.current_step])
            / (((-55. - self.voltage[self.current_step]) / 10.).exp() - 1.)
    }

    fn beta_m(&self) -> f32 {
        4. * ((-self.voltage[self.current_step] + VOLTAGE_INIT) / 18.).exp()
    }

    fn beta_h(&self) -> f32 {
        1. / (((-35. - self.voltage[self.current_step]) / 10.).exp() + 1.)
    }

    fn beta_n(&self) -> f32 {
        0.125 * ((VOLTAGE_INIT - self.voltage[self.current_step]) / 80.).exp()
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

#[cfg(test)]
mod test {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn hodgkin_huxley() {
        let mut hh_model = HodgkinHuxley::new();
        hh_model.time = 1;
        hh_model.dt = 0.001;

        hh_model.calc();

        assert_approx_eq!(hh_model.voltage[1], VOLTAGE_INIT);
        assert_approx_eq!(
            hh_model.i_na[1],
            CONDUCTANCE_NA * GATE_M_INIT.powi(3) * GATE_H_INIT * VOLTAGE_INIT
                - EQUILIBRIUM_VOLTAGE_NA
        );
        assert_approx_eq!(
            hh_model.i_k[1],
            CONDUCTANCE_K * GATE_N_INIT.powi(4) * (VOLTAGE_INIT - EQUILIBRIUM_VOLTAGE_K)
        );
        assert_approx_eq!(
            hh_model.i_leak[1],
            CONDUCTANCE_LEAK * (VOLTAGE_INIT - EQUILIBRIUM_VOLTAGE_LEAK)
        );
        assert_approx_eq!(hh_model.i_injection[1], 30.);
        assert_approx_eq!(
            hh_model.m[1],
            GATE_M_INIT
                + hh_model.diff_gate(
                    0.1 * (-40. - VOLTAGE_INIT) / (((-40. - VOLTAGE_INIT) / 10.).exp() - 1.),
                    4. * 0f32.exp(),
                    GATE_M_INIT
                ) * 0.001
        );
        assert_approx_eq!(
            hh_model.h[1],
            GATE_H_INIT
                + hh_model.diff_gate(
                    0.07 * 0f32.exp(),
                    1. / (((-35. - VOLTAGE_INIT) / 10.).exp() + 1.),
                    GATE_H_INIT
                ) * 0.001
        );
        assert_approx_eq!(
            hh_model.n[1],
            GATE_N_INIT
                + hh_model.diff_gate(
                    0.01 * (-55. - VOLTAGE_INIT) / (((-55. - VOLTAGE_INIT) / 10.).exp() - 1.),
                    0.125 * 0f32.exp(),
                    GATE_N_INIT
                ) * 0.001
        );
    }
}
