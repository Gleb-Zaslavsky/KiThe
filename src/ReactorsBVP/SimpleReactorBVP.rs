use crate::Kinetics::User_reactions::KinData;
use RustedSciThe::numerical::BVP_api::BVP;

use crate::Kinetics::mechfinder_api::kinetics::ElementaryStruct;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;

fn Arr(A: f64, n: f64, E: f64) -> Expr {
    let k = ElementaryStruct::new(vec![A, n, E]);
    return k.K_expr();
}

pub struct SimpleReactorTask {
    pub problem_name: Option<String>,
    pub problem_description: Option<String>,
    pub stoichiometry: KinData,
    pub reaction_data: Option<Vec<HashMap<String, Vec<f64>>>>,
    pub arrhenius_parameters: Vec<(f64, f64, f64)>,
    pub thermal_effects: Vec<f64>,
    pub Cp: f64,
    pub Lambda: f64,
    pub m: f64,
    pub scaling: HashMap<String, f64>,
}

impl SimpleReactorTask {
    pub fn new() -> Self {
        Self {
            problem_name: None,
            problem_description: None,
            stoichiometry: KinData::new(),
            reaction_data: None,
            arrhenius_parameters: Vec::new(),
            thermal_effects: Vec::new(),
            Cp: 0.0,
            Lambda: 0.0,
            m: 0.0,
            scaling: HashMap::new(),
        }
    }

    pub fn set_problem_name(&mut self, name: &str) {
        self.problem_name = Some(name.to_string());
    }

    pub fn set_problem_description(&mut self, description: &str) {
        self.problem_description = Some(description.to_string());
    }
    pub fn set_vector_of_reactions(
        &mut self,
        vector_of_reactions: Vec<String>,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) {
        self.stoichiometry
            .set_reactions_directly(vector_of_reactions, groups);
        self.stoichiometry.analyze_reactions();
    }
    pub fn set_reaction_data(&mut self, reaction_data: Vec<HashMap<String, Vec<f64>>>) {
        self.reaction_data = Some(reaction_data.clone());
        let mut arrhenius_parameters = Vec::new();
        let mut thermal_effects = Vec::new();
        for i in 0..reaction_data.len() {
            let reaction = &reaction_data[i];
            let Q = reaction.get("Q").expect("for reaction Q is not defined")[0];
            thermal_effects.push(Q);
            let kin = reaction
                .get("kin")
                .expect("for reaction Arrhenius parameters is not defined")
                .as_slice();
            arrhenius_parameters.push((kin[0], kin[1], kin[2]));
        }
    }
    pub fn set_arrhenius_parameters(&mut self, arrhenius_parameters: Vec<(f64, f64, f64)>) {
        self.arrhenius_parameters = arrhenius_parameters;
    }
    pub fn set_thermal_effects(&mut self, thermal_effects: Vec<f64>) {
        self.thermal_effects = thermal_effects;
    }
    pub fn parse_from_file(&mut self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(file_path)?;
        let sections: Vec<&str> = content.split("\n\n").collect();

        for section in sections {
            let lines: Vec<&str> = section.lines().collect();
            if lines.is_empty() {
                continue;
            }

            match lines[0].trim() {
                "DESCRIPTION" => {
                    let description = lines[1..].join("\n");
                    self.set_problem_description(&description);
                }
                "REACTIONS" => {
                    let mut reaction_data = Vec::new();
                    let mut vector_of_reactions = Vec::new();
                    for line in &lines[1..] {
                        let parts: Vec<&str> = line.split_whitespace().collect();
                        if parts.len() < 5 {
                            continue;
                        }

                        // Parse the line: reaction Q A n Ea
                        let reaction_eq = parts[0].to_string();
                        let Q: f64 = parts[1].parse()?;
                        let A: f64 = parts[2].parse()?;
                        let n: f64 = parts[3].parse()?;
                        let Ea: f64 = parts[4].parse()?;

                        let mut reaction_map = HashMap::new();

                        reaction_map.insert("Q".to_string(), vec![Q]);
                        reaction_map.insert("kin".to_string(), vec![A, n, Ea]);

                        reaction_data.push(reaction_map);
                        vector_of_reactions.push(reaction_eq)
                    }

                    self.set_reaction_data(reaction_data);
                    self.set_vector_of_reactions(vector_of_reactions, None);
                }
                _ => continue,
            }
        }
        Ok(())
    }
}

pub struct SimpleReactorSolver {
    pub task: SimpleReactorTask,
    pub solver: BVP,
}
