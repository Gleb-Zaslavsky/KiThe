use super::SimpleReactorBVP::SimpleReactorTask;
use RustedSciThe::Utils::logger::{save_matrix_to_csv, save_matrix_to_file};
use RustedSciThe::Utils::plots::{plots, plots_gnulot};
use RustedSciThe::symbolic::symbolic_integration::QuadMethod;
pub struct AfterSolution {}
impl SimpleReactorTask {
    /////////////////////////////////////////POSTPROCESSING///////////////////////////////////////////
    // solver returns scaled variables - now we must return them from dimensionless to dimension form
    pub fn postprocessing(&mut self) {
        if let Some(x_mesh) = self.solver.x_mesh.as_mut() {
            let L = self.L;
            x_mesh.iter_mut().for_each(|xi| *xi *= L);
        }
        if let Some(y) = self.solver.solution.as_mut() {
            let unknowns = &self.solver.unknowns;
            assert_eq!(unknowns.len(), y.ncols());
            let dT = self.scaling.dT;
            let T_scale = self.scaling.T_scale;
            let L = self.L;
            for (i, mut sol_for_var) in y.column_iter_mut().enumerate() {
                if &unknowns[i] == "Teta" {
                    sol_for_var
                        .iter_mut()
                        .for_each(|Teta_i| *Teta_i = *Teta_i * T_scale + dT);
                }
                if &unknowns[i] == "q" {
                    sol_for_var
                        .iter_mut()
                        .for_each(|q_i| *q_i = *q_i * T_scale / L);
                }
                if unknowns[i].starts_with("J") {
                    sol_for_var.iter_mut().for_each(|J_i| *J_i = *J_i / L);
                }
            }
        }
    }
    ////////////////////////////////////////////////I/O/////////////////////////////////////////////////////
    pub fn plot(&self) {
        let y = self.solver.solution.clone().unwrap();
        let x_mesh = self.solver.x_mesh.clone().unwrap();
        let arg = "x".to_owned();
        let values = self.solver.unknowns.clone();
        plots(arg, values, x_mesh, y);
    }
    pub fn gnuplot(&self) {
        let y = self.solver.solution.clone().unwrap();
        let x_mesh = self.solver.x_mesh.clone().unwrap();
        let arg = "x".to_owned();
        let values = self.solver.unknowns.clone();
        plots_gnulot(arg, values, x_mesh, y);
    }

    pub fn save_to_file(&self, filename: Option<String>) {
        let name = if let Some(name) = filename {
            format!("{}.txt", name)
        } else {
            "result.txt".to_string()
        };

        let y = self.solver.solution.clone().unwrap();
        let x_mesh = self.solver.x_mesh.clone().unwrap();
        let _ = save_matrix_to_file(
            &y,
            &self.solver.unknowns,
            &name,
            &x_mesh,
            &self.solver.arg_name,
        );
    }
    pub fn save_to_csv(&self, filename: Option<String>) {
        let name = if let Some(name) = filename {
            name
        } else {
            "result_table".to_string()
        };
        let y = self.solver.solution.clone().unwrap();
        let x_mesh = self.solver.x_mesh.clone().unwrap();
        let _ = save_matrix_to_csv(
            &y,
            &self.solver.unknowns,
            &name,
            &x_mesh,
            &self.solver.arg_name,
        );
    }

    pub fn estimate_values(&self) {
        if self.kindata.vec_of_equations.len() == 1 {
            let Q = self.thermal_effects[0];
            let C_p = self.Cp;
            let T0 = self
                .boundary_condition
                .get("T")
                .expect("No T in boundary condition");
            let T_fin = T0 + Q / C_p;
            println!(
                " if there will be an adiabatic reactor of ideal mixing the temperature
            will be equal to {}",
                T_fin
            );
        }
        let unknowns = self.solver.unknowns.clone();
        let heat_release = self.heat_release.clone();
        let x = self.solver.x_mesh.clone().unwrap();

        let y = self.solver.solution.clone().unwrap();
        let heat_release_fun = heat_release.lambdify(unknowns.iter().map(|x| x.as_str()).collect());
    }

    fn energy_balance(&self) {
        let i = self.solver.unknowns.iter().position(|x| x == "q").unwrap();
        let q:Vec<f64> = self.solver.solution.as_ref().unwrap().column(i).iter().cloned().collect::<Vec<f64>>();
        let q_f = q.last().unwrap();
        let q_0 = q.first().unwrap();
        let dq = q_f - q_0;
        let i = self.solver.unknowns.iter().position(|x| x == "Teta").unwrap();
        let T = self.solver.solution.as_ref().unwrap().column(i).iter().cloned().collect::<Vec<f64>>();
        let T_f = T.last().unwrap();
        let T_0 = T.first().unwrap();
        let m = self.m;
        let Cp = self.Cp;
        let dT = m*Cp*(T_f - T_0);


    }
}
