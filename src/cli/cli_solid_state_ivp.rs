use crate::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
use RustedSciThe::numerical::ODE_api2::SolverType;
use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
use dialoguer::{Input, Select, theme::ColorfulTheme};

pub fn solid_state_ivp_menu() {
    let theme = ColorfulTheme::default();
    loop {
        let items = vec!["Solve kinetic problem", "Exit"];
        let selection = Select::with_theme(&theme)
            .with_prompt("=== Solid State Kinetics IVP Solver ===")
            .items(&items)
            .default(0)
            .interact()
            .unwrap();

        match selection {
            0 => {
                if let Err(e) = run_solver() {
                    println!("Error: {}", e);
                }
            }
            1 => break,
            _ => {}
        }
    }
}

fn run_solver() -> Result<(), String> {
    // Step 1: Choose solver type
    let solver_type = choose_solver_type()?;

    // Step 2: Create solver instance
    let mut ivp = KineticModelIVP::new(solver_type);

    // Step 3: Choose kinetic model
    let (model_name, params) = choose_kinetic_model()?;

    // Step 4: Input parameters
    let (t_final, beta, t0, e, a) = input_parameters()?;

    // Step 5: Set up and solve
    ivp.set_problem(t_final, beta, t0, e, a)?;
    ivp.set_model(model_name, params)?;

    println!("Solving...");
    ivp.solve()?;

    println!("Solution complete! Plotting results...");
    ivp.plot()?;
    ivp.plot_in_terminal();
    Ok(())
}

fn choose_solver_type() -> Result<SolverType, String> {
    let theme = ColorfulTheme::default();
    let items = vec![
        "BDF (recommended for stiff problems!!!OPTION NON PRODUCT READY)",
        "Radau (high accuracy)",
        "RK45 (non-stiff problems)",
        "Backward Euler (simple implicit)",
    ];
    let selection = Select::with_theme(&theme)
        .with_prompt("Choose solver type")
        .items(&items)
        .default(1)
        .interact()
        .unwrap();

    match selection {
        0 => Ok(SolverType::BDF),
        1 => Ok(SolverType::Radau(RadauOrder::Order3)),
        2 => Ok(SolverType::NonStiff("RK45".to_owned())),
        3 => Ok(SolverType::BackwardEuler),
        _ => Err("Invalid solver choice".to_string()),
    }
}

fn choose_kinetic_model() -> Result<(KineticModelNames, Vec<f64>), String> {
    println!("\nAvailable kinetic models:");
    KineticModelNames::pretty_print();

    let theme = ColorfulTheme::default();
    let model_code: String = Input::with_theme(&theme)
        .with_prompt("Enter model code (e.g., A2, JMA, SB)")
        .interact_text()
        .unwrap();
    let model_code = model_code.trim().to_uppercase();

    let model_name = match model_code.as_str() {
        "A2" => KineticModelNames::A2,
        "A3" => KineticModelNames::A3,
        "A4" => KineticModelNames::A4,
        "R2" => KineticModelNames::R2,
        "R3" => KineticModelNames::R3,
        "D1" => KineticModelNames::D1,
        "D2" => KineticModelNames::D2,
        "D3" => KineticModelNames::D3,
        "D4" => KineticModelNames::D4,
        "P2_3" => KineticModelNames::P2_3,
        "P2" => KineticModelNames::P2,
        "P3" => KineticModelNames::P3,
        "F1" => KineticModelNames::F1,
        "F2" => KineticModelNames::F2,
        "F3" => KineticModelNames::F3,
        "SB" => KineticModelNames::SB,
        "JMA" => KineticModelNames::JMA,
        "AC" => KineticModelNames::Ac,
        "DEC" => KineticModelNames::Dec,
        "PTE" => KineticModelNames::PTe,
        "SBTP" => KineticModelNames::SBtp,
        _ => return Err("Invalid model code".to_string()),
    };

    let params = match model_name {
        KineticModelNames::SB => {
            println!("SB model requires 3 parameters (m, n, p):");
            input_model_params(3)?
        }
        KineticModelNames::JMA => {
            println!("JMA model requires 1 parameter (m):");
            input_model_params(1)?
        }
        KineticModelNames::Ac => {
            println!("Ac model requires 1 parameter (m):");
            input_model_params(1)?
        }
        KineticModelNames::Dec => {
            println!("Dec model requires 1 parameter (m):");
            input_model_params(1)?
        }
        KineticModelNames::PTe => {
            println!("PTe model requires 2 parameters (m, n):");
            input_model_params(2)?
        }
        KineticModelNames::SBtp => {
            println!("SBtp model requires 3 parameters (m, n, c):");
            input_model_params(3)?
        }
        _ => vec![], // Parameter-free models
    };

    Ok((model_name, params))
}

fn input_model_params(count: usize) -> Result<Vec<f64>, String> {
    let theme = ColorfulTheme::default();
    let mut params = Vec::new();

    for i in 0..count {
        let prompt = format!("Parameter {}", i + 1);
        let val: f64 = Input::<f64>::with_theme(&theme)
            .with_prompt(prompt)
            .interact_text()
            .map_err(|_| "Invalid number format".to_string())?;
        params.push(val);
    }

    Ok(params)
}

fn input_parameters() -> Result<(f64, f64, f64, f64, f64), String> {
    println!("\nEnter problem parameters:");
    let theme = ColorfulTheme::default();

    let t_final: f64 = Input::<f64>::with_theme(&theme)
        .with_prompt("t_final (final time, s)")
        .interact_text()
        .map_err(|_| "Invalid t_final format".to_string())?;

    let beta: f64 = Input::<f64>::with_theme(&theme)
        .with_prompt("beta (heating rate, K/s)")
        .interact_text()
        .map_err(|_| "Invalid beta format".to_string())?;

    let t0: f64 = Input::<f64>::with_theme(&theme)
        .with_prompt("T0 (initial temperature, K)")
        .interact_text()
        .map_err(|_| "Invalid T0 format".to_string())?;

    let e: f64 = Input::<f64>::with_theme(&theme)
        .with_prompt("E (activation energy, J/mol)")
        .interact_text()
        .map_err(|_| "Invalid E format".to_string())?;

    let a: f64 = Input::<f64>::with_theme(&theme)
        .with_prompt("A (pre-exponential factor, 1/s)")
        .interact_text()
        .map_err(|_| "Invalid A format".to_string())?;

    Ok((t_final, beta, t0, e, a))
}
