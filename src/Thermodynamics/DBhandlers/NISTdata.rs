use std::collections::HashMap;
use std::f64::consts::LN_2;

fn log(x: f64) -> f64 {
    x.ln()
}

fn calculate_cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    a + b * t + c * t.powi(2) + d * t.powi(3) + e / t.powi(2)
}

fn calculate_dh(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, g: f64, h: f64) -> f64 {
    a * t + (b * t.powi(2)) / 2.0 + (c * t.powi(3)) / 3.0 + (d * t.powi(4)) / 4.0 - e / t + f - h
}

fn calculate_s(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, g: f64, h: f64) -> f64 {
    a * log(t) + b * t + (c * t.powi(2)) / 2.0 + (d * t.powi(3)) / 3.0 - e / (2.0 * t.powi(2)) + g
}
/*
fn preliminary_work_with_cp(cp_data: &Vec<f64>) -> Vec<f64> {
    let condition = if cp_data[0..2] == cp_data[2..4] && cp_data[2..4] == cp_data[4..6] {
        11
    } else if cp_data[0..2] == cp_data[2..4] && cp_data[2..4] != cp_data[4..6] {
        12
    } else if cp_data[0..2] != cp_data[2..4] && cp_data[2..4] == cp_data[4..6] {
        21
    } else {
        22
    };

    let (cp_data_start, cp_data_end) = match condition {
        11 => (cp_data[0..2].to_vec(), cp_data[6..].to_vec()),
        12 => (vec![cp_data[0], cp_data[1], cp_data[4]], cp_data[5..].to_vec()),
        21 => {
            if cp_data[1] == cp_data[2] {
                (vec![cp_data[0], cp_data[1], cp_data[3]], cp_data[6..].to_vec())
            } else {
                (cp_data[0..4].to_vec(), cp_data[6..].to_vec())
            }
        }
        22 => {
            if cp_data[1] == cp_data[2] && cp_data[3] == cp_data[4] {
                (vec![cp_data[0], cp_data[1], cp_data[3]], cp_data[5..].to_vec())
            } else if cp_data[1] == cp_data[2] && cp_data[3] != cp_data[4] {
                (vec![cp_data[0], cp_data[1], cp_data[3]], cp_data[4..].to_vec())
            } else {
                (cp_data[0..6].to_vec(), cp_data[6..].to_vec())
            }
        }
        _ => (vec![], vec![]),
    };

    [cp_data_start, cp_data_end].concat()
}

fn calc_constants_nist(
    tm: f64,
    data_substance: &HashMap<&str, f64>,
    cp_from_nist: Vec<f64>,
    flag_calc: i32,
) -> HashMap<&'static str, f64> {
    let t = tm / 1000.0;
    let mut nist_constant: HashMap<&str, f64> = HashMap::new();
    nist_constant.insert("Cp", 0.0);
    nist_constant.insert("dH", 0.0);
    nist_constant.insert("dS", 0.0);

    let mut cp_data = if flag_calc == 1 {
        data_substance.get("Cp").unwrap().to_vec()
    } else {
        cp_from_nist
    };

    if cp_data.len() != 1 {
        cp_data = preliminary_work_with_cp(&cp_data);
    }

    if cp_data.len() != 1 {
        if cp_data.len() == 28 {
            let t1 = cp_data[0];
            let t2 = cp_data[1];
            let t3 = cp_data[2];
            let t4 = cp_data[3];
            if t1 <= tm && tm <= t2 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[4], cp_data[5], cp_data[6], cp_data[7], cp_data[8], cp_data[9], cp_data[10], cp_data[11],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else if t2 < tm && tm <= t3 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[12], cp_data[13], cp_data[14], cp_data[15], cp_data[16], cp_data[17], cp_data[18], cp_data[19],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else if t3 < tm && tm < t4 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[20], cp_data[21], cp_data[22], cp_data[23], cp_data[24], cp_data[25], cp_data[26], cp_data[27],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else {
                nist_constant.insert("Cp", f64::NAN);
            }
        } else if cp_data.len() == 19 {
            let t1 = cp_data[0];
            let t2 = cp_data[1];
            let t3 = cp_data[2];
            if t1 <= tm && tm <= t2 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[3], cp_data[4], cp_data[5], cp_data[6], cp_data[7], cp_data[8], cp_data[9], cp_data[10],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else if t2 < tm && tm <= t3 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[11], cp_data[12], cp_data[13], cp_data[14], cp_data[15], cp_data[16], cp_data[17], cp_data[18],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else {
                nist_constant.insert("Cp", f64::NAN);
            }
        } else if cp_data.len() == 10 {
            let t1 = cp_data[0];
            let t2 = cp_data[1];
            if t1 <= tm && tm <= t2 {
                let (a, b, c, d, e, f, g, h) = (
                    cp_data[2], cp_data[3], cp_data[4], cp_data[5], cp_data[6], cp_data[7], cp_data[8], cp_data[9],
                );
                if flag_calc == 1 {
                    let cp = calculate_cp(t, a, b, c, d, e);
                    nist_constant.insert("Cp", cp);
                } else {
                    let dh_0 = data_substance["dH"];
                    let dh1 = calculate_dh(t, a, b, c, d, e, f, g, h);
                    let dh = 1000.0 * (dh_0 + dh1);
                    let ds = calculate_s(t, a, b, c, d, e, f, g, h);
                    nist_constant.insert("dH", dh);
                    nist_constant.insert("dS", ds);
                }
            } else {
                nist_constant.insert("Cp", f64::NAN);
            }
        }
    } else if cp_data.len() == 1 {
        let cp = cp_data[0];
        let dh_0 = data_substance["dH"];
        let dh1 = cp * (tm - 298.0);
        let ds = data_substance["dS"];
        let dh = 1000.0 * dh_0 + dh1;
        nist_constant.insert("Cp", cp);
        nist_constant.insert("dH", dh);
        nist_constant.insert("dS", ds);
    } else {
        nist_constant.insert("Cp", 0.0);
        nist_constant.insert("dH", 0.0);
        nist_constant.insert("dS", 0.0);
    }

    nist_constant
}

fn ma() {
    let tm = 300.0;
    let mut data_substance = HashMap::new();
    data_substance.insert("Cp", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    data_substance.insert("dH", 100.0);
    data_substance.insert("dS", 200.0);
    let cp_from_nist = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let flag_calc = 1;

    let nist_constant = calc_constants_nist(tm, &data_substance, cp_from_nist, flag_calc);
    println!("{:?}", nist_constant);
}
 */
