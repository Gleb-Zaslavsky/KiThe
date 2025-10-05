use crate::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
use RustedSciThe::numerical::ODE_api2::SolverType;
use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
use std::io::{self, Write};

pub fn solid_state_ivp_menu() {
    loop {
        println!("\n=== Solid State Kinetics IVP Solver ===");
        println!("1. Solve kinetic problem");
        println!("2. Exit");
        print!("Choose option: ");
        io::stdout().flush().unwrap();

        let mut input = String::new();
        io::stdin().read_line(&mut input).unwrap();

        match input.trim() {
            "1" => {
                if let Err(e) = run_solver() {
                    println!("Error: {}", e);
                }
            }
            "2" => break,
            _ => println!("Invalid option"),
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
    println!("\nChoose solver type:");
    println!("1. BDF (recommended for stiff problems!!!OPTION NON PRODUCT READY)");
    println!("2. Radau (high accuracy)");
    println!("3. RK45 (non-stiff problems)");
    println!("4. Backward Euler (simple implicit)");
    print!("Enter choice (1-4): ");
    io::stdout().flush().unwrap();

    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();

    match input.trim() {
        "1" => Ok(SolverType::BDF),
        "2" => Ok(SolverType::Radau(RadauOrder::Order3)),
        "3" => Ok(SolverType::NonStiff("RK45".to_owned())),
        "4" => Ok(SolverType::BackwardEuler),
        _ => Err("Invalid solver choice".to_string()),
    }
}

fn choose_kinetic_model() -> Result<(KineticModelNames, Vec<f64>), String> {
    println!("\nAvailable kinetic models:");
    KineticModelNames::pretty_print();

    print!("\nEnter model code (e.g., A2, JMA, SB): ");
    io::stdout().flush().unwrap();

    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    let model_code = input.trim().to_uppercase();

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
    let mut params = Vec::new();

    for i in 0..count {
        print!("Parameter {}: ", i + 1);
        io::stdout().flush().unwrap();

        let mut input = String::new();
        io::stdin().read_line(&mut input).unwrap();

        let param: f64 = input
            .trim()
            .parse()
            .map_err(|_| "Invalid number format".to_string())?;
        params.push(param);
    }

    Ok(params)
}

fn input_parameters() -> Result<(f64, f64, f64, f64, f64), String> {
    println!("\nEnter problem parameters:");

    print!("t_final (final time, s): ");
    io::stdout().flush().unwrap();
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    let t_final: f64 = input
        .trim()
        .parse()
        .map_err(|_| "Invalid t_final format".to_string())?;

    print!("beta (heating rate, K/s): ");
    io::stdout().flush().unwrap();
    input.clear();
    io::stdin().read_line(&mut input).unwrap();
    let beta: f64 = input
        .trim()
        .parse()
        .map_err(|_| "Invalid beta format".to_string())?;

    print!("T0 (initial temperature, K): ");
    io::stdout().flush().unwrap();
    input.clear();
    io::stdin().read_line(&mut input).unwrap();
    let t0: f64 = input
        .trim()
        .parse()
        .map_err(|_| "Invalid T0 format".to_string())?;

    print!("E (activation energy, J/mol): ");
    io::stdout().flush().unwrap();
    input.clear();
    io::stdin().read_line(&mut input).unwrap();
    let e: f64 = input
        .trim()
        .parse()
        .map_err(|_| "Invalid E format".to_string())?;

    print!("A (pre-exponential factor, 1/s): ");
    io::stdout().flush().unwrap();
    input.clear();
    io::stdin().read_line(&mut input).unwrap();
    let a: f64 = input
        .trim()
        .parse()
        .map_err(|_| "Invalid A format".to_string())?;

    Ok((t_final, beta, t0, e, a))
}
