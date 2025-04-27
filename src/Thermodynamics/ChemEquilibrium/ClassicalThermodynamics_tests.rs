#[cfg(test)]
mod tests {

    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};

    use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::{
        ThermodynamicCalculations, Thermodynamics,
    };

    use crate::Thermodynamics::User_PhaseOrSolution::{
        CustomSubstance, PhaseOrSolution, SubstanceSystemFactory, SubstancesContainer,
    };
    use approx::assert_relative_eq;
    use std::collections::HashMap;

    #[test]
    fn test_substances_container() {
        // Test SinglePhase container
        let substances = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(substances.clone());

        let all_substances = container.get_all_substances();
        assert_eq!(all_substances.len(), 2);
        assert!(all_substances.contains(&"CO".to_string()));
        assert!(all_substances.contains(&"CO2".to_string()));

        // Test MultiPhase container
        let mut phase_substances = HashMap::new();
        phase_substances.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phase_substances.insert("liquid".to_string(), vec!["H2O".to_string()]);

        let container = SubstancesContainer::MultiPhase(phase_substances);
        let all_substances = container.get_all_substances();

        assert_eq!(all_substances.len(), 3);
        assert!(all_substances.contains(&"H2".to_string()));
        assert!(all_substances.contains(&"O2".to_string()));
        assert!(all_substances.contains(&"H2O".to_string()));
    }

    #[test]
    fn test_phase_or_solution_new() {
        let phase_or_solution = PhaseOrSolution::new();
        assert!(phase_or_solution.subs_data.is_empty());
    }

    #[test]
    fn test_custom_substance_one_phase() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string(), "CO2".to_string()];

        // Set library priorities
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Set phases
        subs_data
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        subs_data
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));

        let custom_substance = CustomSubstance::OnePhase(subs_data);

        // Test that the custom substance was created correctly
        match &custom_substance {
            CustomSubstance::OnePhase(data) => {
                assert_eq!(data.substances.len(), 2);
                assert!(data.substances.contains(&"CO".to_string()));
                assert!(data.substances.contains(&"CO2".to_string()));
            }
            _ => panic!("Expected OnePhase variant"),
        }
    }

    #[test]
    fn test_substance_system_factory_single_phase() {
        let substances = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        assert!(result.is_ok());

        let custom_substance = result.unwrap();
        match custom_substance {
            CustomSubstance::OnePhase(data) => {
                assert_eq!(data.substances.len(), 2);
                assert!(data.map_of_phases.contains_key("CO"));
                assert!(data.map_of_phases.contains_key("CO2"));
            }
            _ => panic!("Expected OnePhase variant"),
        }
    }

    #[test]
    fn test_substance_system_factory_multi_phase() {
        let mut phase_substances = HashMap::new();
        phase_substances.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phase_substances.insert("liquid".to_string(), vec!["H2O".to_string()]);

        let container = SubstancesContainer::MultiPhase(phase_substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        assert!(result.is_ok());

        let custom_substance = result.unwrap();
        match custom_substance {
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                assert_eq!(phase_or_solution.subs_data.len(), 2);
                assert!(
                    phase_or_solution
                        .subs_data
                        .contains_key(&Some("gas".to_string()))
                );
                assert!(
                    phase_or_solution
                        .subs_data
                        .contains_key(&Some("liquid".to_string()))
                );

                // Check gas phase substances
                let gas_data = phase_or_solution
                    .subs_data
                    .get(&Some("gas".to_string()))
                    .unwrap();
                assert_eq!(gas_data.substances.len(), 2);
                assert!(gas_data.substances.contains(&"H2".to_string()));
                assert!(gas_data.substances.contains(&"O2".to_string()));

                // Check liquid phase substances
                let liquid_data = phase_or_solution
                    .subs_data
                    .get(&Some("liquid".to_string()))
                    .unwrap();
                assert_eq!(liquid_data.substances.len(), 1);
                assert!(liquid_data.substances.contains(&"H2O".to_string()));
            }
            _ => panic!("Expected PhaseOrSolution variant"),
        }
    }

    #[test]
    fn test_thermodynamic_calculations_trait() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string(), "CO2".to_string()];

        // Set library priorities
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Set phases
        subs_data
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        subs_data
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));

        // Search for substance data
        subs_data.search_substances();

        let mut custom_substance = CustomSubstance::OnePhase(subs_data);

        // Test create_thermodynamics method
        let mut n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)> = HashMap::new();
        n.insert(None, (Some(1.0), Some(vec![0.5, 0.5])));
        let result = custom_substance.create_thermodynamics(400.0, 101325.0, Some(n), None);

        assert!(result.is_ok());

        let thermodynamics = result.unwrap();
        assert_eq!(thermodynamics.T, 400.0);
        assert_eq!(thermodynamics.P, 101325.0);
        assert!(!thermodynamics.dG.is_empty());
        assert!(!thermodynamics.dG_sym.is_empty());
        assert!(!thermodynamics.dG_fun.is_empty());
    }

    #[test]
    fn test_thermodynamics_new() {
        let thermo = Thermodynamics::new();
        assert_eq!(thermo.T, 298.15);
        assert_eq!(thermo.P, 1e5);
        assert!(thermo.dG.is_empty());
        assert!(thermo.dG_sym.is_empty());
        assert!(thermo.dG_fun.is_empty());
    }

    #[test]
    fn test_thermodynamics_setters() {
        let mut thermo = Thermodynamics::new();

        // Test set_T
        thermo.set_T(400.0);
        assert_eq!(thermo.T, 400.0);

        // Test set_P
        thermo.set_P(200000.0, None);
        assert_eq!(thermo.P, 200000.0);
    }

    #[test]
    fn test_thermodynamics_update_substances() {
        // Create a test SubsData
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string(), "CO2".to_string()];

        // Create a Thermodynamics instance with the test SubsData
        let mut thermo = Thermodynamics::new();
        thermo.subdata = CustomSubstance::OnePhase(subs_data);

        // Update substances from subdata
        thermo.update_substances_from_subdata();

        // Check that substances were updated correctly
        match &thermo.subs_container {
            SubstancesContainer::SinglePhase(substances) => {
                assert_eq!(substances.len(), 2);
                assert!(substances.contains(&"CO".to_string()));
                assert!(substances.contains(&"CO2".to_string()));
            }
            _ => panic!("Expected SinglePhase variant"),
        }
    }

    #[test]
    fn test_thermodynamics_calculations_integration() {
        // Create a test system using the factory
        let substances = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        assert!(result.is_ok());

        let mut custom_substance = result.unwrap();
        let mut n = HashMap::new();
        n.insert(None, (Some(1.0), Some(vec![0.5, 0.5])));
        // Create thermodynamics from the custom substance
        let result = custom_substance.create_thermodynamics(400.0, 101325.0, Some(n.clone()), None);

        assert!(result.is_ok());

        let mut thermo = result.unwrap();

        // Test that calculations were performed correctly
        assert!(!thermo.dG.is_empty());
        assert!(!thermo.dG_sym.is_empty());
        assert!(!thermo.dS.is_empty());
        assert!(!thermo.dS_sym.is_empty());

        // Test recalculation at a different temperature
        thermo.calculate_Gibbs_free_energy(500.0, n);
        assert_eq!(thermo.T, 500.0);

        // Test set_P_to_sym
        thermo.set_P_to_sym();
        // The pressure should be substituted in the symbolic expressions
    }

    #[test]
    fn test_thermodynamics_clone() {
        // Create a test system
        let substances = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        let mut custom_substance = result.unwrap();

        // Create thermodynamics
        let result = custom_substance.create_thermodynamics(400.0, 101325.0, None, None);

        let thermo = result.unwrap();

        // Clone the thermodynamics instance
        let thermo_clone = thermo.clone();

        // Test that the clone has the same values
        assert_eq!(thermo_clone.T, thermo.T);
        assert_eq!(thermo_clone.P, thermo.P);

        // Check that dG maps are the same
        for (phase, substances) in &thermo.dG {
            assert!(thermo_clone.dG.contains_key(phase));
            let clone_substances = thermo_clone.dG.get(phase).unwrap();

            for (substance, value) in substances {
                assert!(clone_substances.contains_key(substance));
                assert_relative_eq!(
                    clone_substances.get(substance).unwrap(),
                    value,
                    epsilon = 1e-6
                );
            }
        }

        // Similarly for dG_sym
        for (phase, substances) in &thermo.dG_sym {
            assert!(thermo_clone.dG_sym.contains_key(phase));
            let clone_substances = thermo_clone.dG_sym.get(phase).unwrap();

            for (substance, _) in substances {
                assert!(clone_substances.contains_key(substance));
                // We can't directly compare Expr objects, but we can check they exist
            }
        }
    }
}
