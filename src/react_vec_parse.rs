 
/// ru
/// Модуль берет на вход  вектор уравнений реации, заданных в виде String и выдает следующие данные: 
/// 1) стехиометрическую матрицу заданную в виде вектора векторов 
/// 2) вектор веществ 
/// 3) вектор векторов стехиометрических коэффициентов реагентов в каждой реакции 
/// 4) то же для продуктов
/// В процессе производится избавление от артефактов парсинга уравнений из баз данных
/// ----------------------------------------------------------------
/// eng
/// The module takes as input a vector of reaction equations specified as a vector of String and produces the following data:
 /// 1) a stoichiometric matrix specified as a vector of vectors 
 /// 2) a vector of substances 
 /// 3) a vector of vectors of stoichiometric coefficients of reactants in each reaction 
 /// 4) the same for products
/// The process involves getting rid of artifacts of parsing equations from databases
/// # Examples
/// reaction data with  artifacts of parsing equations from databases
/// ```
/// let reactions_: Vec<&str> = vec!["A=2BM)", "B->A + 3C_DUP", "2B+A=D"];
/// let reaction = reactions_.iter().map(|s| s.to_string()).collect();
/// ReactionAnalyzer_instance.reactions = reaction;
/// ReactionAnalyzer_instance.search_substances();
/// println!("substances: {:?}", ReactionAnalyzer_instance.substances);
/// let substancses_ = vec!["A", "B", "C", "D"];
/// let substancses = substancses_.iter().map(|s| s.to_string()).collect();
 /// ReactionAnalyzer_instance.substances = substancses ;
/// ReactionAnalyzer_instance.analyse_reactions();
/// println!("{:?}", ReactionAnalyzer_instance);
/// ```
use regex::Regex;
use std::collections::HashSet;
#[warn(unused_variables)]
#[warn(unused_mut)]



/*
        #  строковая переменная на которые разделяет по знаку "+" половину уравнения относящююся к 
        # реагентам/продуктам может содержать '(', ')', а также 'M' - символ третьей частицы которая не относится к 
        # продуктам и реагентам
*/
fn clean_off_artifacts( item: &mut String) -> &mut String {

    let items_to_clean = [" (", "(", "M)", "M)"];
    if item.contains( "M)") || item.contains("M)") || item.contains(" (") || item.contains("(") {
        *item = item.replace("M)", "").replace("M)", "").replace(" (", "").replace("(", "");
        
        return item
    } else {item}
}

pub fn analyse_substances(half_reaction: &str) -> (Vec<String>, Vec<f64>, Vec<f64>) {
    let mut s_list = Vec::new();
    let mut g_list = Vec::new();
    let mut subs:Vec<String> = Vec::new();

    let re_coeff = Regex::new(r"(\d+(\.\d*)?)\*?").unwrap();
    let re_power = Regex::new(r"\^\ *(\d+(\.\d*)?)").unwrap();

    for s in half_reaction.split('+') {
        let mut s = s.trim();
        print!("half reaction is:  {:#?}   ",s);
        let mut stec_coeff = 1.0;
        let mut power_coeff = 1.0;
        let end_of_stoichiometric: usize = 0; 
        //  captures.get(1) is used to access the first capture group in the regular expression match. 
        if let Some(captures) = re_coeff.captures(s) {
            
            // refers to the first capture group in the regular expression match, example: Match { start: 0, end: 1, string: "2" }
            let regmatch = captures.get(0).unwrap();
            print!(" regmatch:  {:?}   ", regmatch);
            let start_of_stoichiometric = regmatch.start();
            let end_of_stoichiometric = regmatch.end();
            // стехиометрический коэффициент должен быть перед формулой вещества, т.е. номер его позиции должен начинаться с нуля
            if start_of_stoichiometric==0{
                s = &s[end_of_stoichiometric..]; 
                stec_coeff = captures.get(1).unwrap().as_str().parse().unwrap();
                print!(" stec_coeff:  {}   ", &stec_coeff);
            }
            
        
        }
    
        

        if let Some(captures) = re_power.captures(s) {
            power_coeff = captures.get(1).unwrap().as_str().parse().unwrap();
            let end = captures.get(0).unwrap().end();
      
        }
        else{
            power_coeff==stec_coeff;
        }



        s = s.trim();
        let s_to_filter: &mut String = &mut s.to_string();
        let s_filtered =  clean_off_artifacts(s_to_filter);
        let s: String = s_filtered.to_string();
        s_list.push(stec_coeff);
        g_list.push(power_coeff);
        subs.push(s);
    }

    (subs, s_list, g_list)
}

// #  некоторые записи уравнения реакции содержат последние символы '_dup' или '_DUP' (то есть дублирующие)
// избавляемся от них 
fn clean_off_DUP( item: &mut String) -> &String {

  
    if item.ends_with("_dup") || item.ends_with("_DUP") {
        *item = item.replace("_dup", "").replace("_DUP", "");
        return item
    } else {item}
    
}

// This struct is used to represent a reaction.
#[derive(Debug)]
pub struct ReactionAnalyzer {
  
    pub reactions: Vec<String>, // вектор реакций/ a vector of reactions
    pub substances: Vec<String>, // вектор веществ/ a vector of substances
    pub stecheo_matrx: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов в каждой реакции/ a vector of vectors of stoichiometric coefficients in each reaction
    pub stecheo_reags: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов реагентов в каждой реакции/ a vector of vectors of stoichiometric coefficients of reactants in each reaction  
    pub stecheo_prods: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов продуктов в каждой реакции/ a vector of vectors of stoichiometric coefficients of products in each reaction
    pub G_matrix_reag: Vec<Vec<f64>>,
    pub G_matrix_prod: Vec<Vec<f64>>,
}

impl ReactionAnalyzer {
    pub fn new() -> Self {
        ReactionAnalyzer {
       
            reactions: Vec::new(),
            substances: Vec::new(),
            stecheo_matrx: Vec::new(),
            stecheo_reags: Vec::new(),
            stecheo_prods: Vec::new(),
            G_matrix_reag: Vec::new(),
            G_matrix_prod: Vec::new(),
        }
    }

    pub fn analyse_reactions(&mut self) {
        let reactions_trimmed = self.reactions.iter().map(|s| s.replace("**", "^").trim().to_string()).collect();
        self.reactions = reactions_trimmed;
        self.stecheo_matrx.resize(self.reactions.len(), vec![0.0; self.substances.len()]);
        self.stecheo_reags.resize(self.reactions.len(), vec![0.0; self.substances.len()]);
        self.stecheo_prods.resize(self.reactions.len(), vec![0.0; self.substances.len()]);
        self.G_matrix_reag.resize(self.reactions.len(), vec![0.0; self.substances.len()]);
        self.G_matrix_prod.resize(self.reactions.len(), vec![0.0; self.substances.len()]);
        for mut reaction in &self.reactions {
            // номер реакции
            if let Some(i) = self.reactions.iter().position(|s| s == reaction){
            println!("reaction number: {}", i);
            let mut reaction_: &mut String = &mut reaction.to_string();
            reaction = clean_off_DUP( reaction_);
            println!("Reaction after dup: {}", &reaction);
            // разделяем уравнение реакции на половины относящиеся к реагентам и продуктам по соответствующему знаку = или -> или =>
            let re = Regex::new(r"=|->|=>").unwrap();
            let sides: Vec<String> = re.split(reaction)// reaction.split(|s| s == '=' ) 
            .map(|s| s.trim())
            .map(|s| s.to_string())
            .collect();
            let mut subs = Vec::new();
            let mut subs_r = Vec::new();
            let mut s_list = Vec::new();
            let mut g_list = Vec::new();
            let mut s_list_r = Vec::new();
            let mut g_list_r: Vec<f64> = Vec::new();
        
            println!("direct reaction  {:?}", &sides[0]);
            let (mut left_subs, mut left_s_list, mut left_g_list) = analyse_substances( &sides[0]);
            let mut left_subs = left_subs.iter().map(|s| s.as_str()).collect();
            subs.append(&mut left_subs);
            s_list.append(&mut left_s_list);
            g_list.append(&mut left_g_list);
            println!("direct reaction substances: {:?}", subs);
            println!("iterating over substances in direct reaction");
            println!("lengh of substance list found in reaction {}", &subs.len());
            for j in 0..subs.len() {
                let subs_j = subs[j];
                if let Some(k) = self.substances.iter().position(|s| s == subs_j) {
                    println!("Index of substance '{}' in list of substanced of this react is: {}, in general list is {},
                     lengh of gen. list is {}", subs_j, j, k,  &self.substances.len());
                    self.stecheo_matrx[i][k] -= s_list[j];
                    self.stecheo_reags[i][k] = s_list[j];
                    self.G_matrix_reag[i][k] = g_list[j];
   
                } else {
                    println!("'{}' not found in the vector", subs_j);
                }



            }
           
            println!("reverse reaction {:?}", &sides[1]);
            let (mut right_subs, mut right_s_list, mut right_g_list) = analyse_substances( &sides[1]);
            let mut right_subs = right_subs.iter().map(|s| s.as_str()).collect();
            subs_r.append(&mut right_subs);
            s_list_r.append(&mut right_s_list);
            g_list_r.append(&mut right_g_list);
            println!("reverse reaction substances: {:?}", subs_r);
            println!("iterating over substances in reverse reaction");
            println!("lengh of substance list found in reaction {}", &subs.len());
            for j in 0..subs_r.len() {          
                let subs_j = subs_r[j];
                // индекс реагента в векторе реагентов
                if let Some(k) = self.substances.iter().position(|s| s == subs_j) {
                    println!("Index of substance '{}' in list of substanced of this react is: {}, in general list is {},
                    lengh of gen. list is {}", subs_j, j, k,  &self.substances.len()); 
                    self.stecheo_matrx[i][k] += s_list_r[j];
                    self.stecheo_prods[i][k] = s_list_r[j];
                    self.G_matrix_prod[i][k] = g_list_r[j];
             
                } else {
                    println!("'{}' not found in the vector", subs_j);
                }



            }
            // Add code to populate self.stecheo_matrx, self.G_matrix_reag, self.G_matrix_prod,
            // self.stecheo_reags, and self.stecheo_prods with the appropriate values.
        } 
   
        } // end of for reaction in &self.reactions {
    }

    pub fn search_substances(&mut self) {
        let reactions_trimmed = self.reactions.iter().map(|s| s.replace("**", "^").trim().to_string()).collect();
        self.reactions = reactions_trimmed;
        let mut found_substances:Vec<String> = Vec::new();
        for mut reaction in &mut self.reactions {
            println!("Reaction: {}", reaction);
          //  let mut reaction: &mut String = &mut reaction.to_string();
            let reaction = clean_off_DUP(reaction);
            let re = Regex::new(r"=|->|=>").unwrap();
            let sides: Vec<&str> = re.split(reaction).map(|s| s).collect();  
            println!("Sides: {:?}", &sides);
            for side in sides {
              let (subs, _, _) = analyse_substances( &side); 
              let subs_mut = &mut subs.clone();
              subs_mut.retain(|s|!found_substances.contains(s));  // remove duplicates
              found_substances.extend(subs_mut.to_owned());
        }
        self.substances = found_substances.iter().map(|s| s.to_string()).collect();
        println!("Substances found: {:?}", &self.substances);

    }
}
}   // end of impl ReactionAnalyzer


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analyse_substances() {
        let half_reaction = "5H2O + 10O2";
        let (subs, s_list, g_list) = analyse_substances(half_reaction);
        assert_eq!(subs, vec!["5H2O".to_string(), "10O2".to_string()]);
        assert_eq!(s_list, vec![5.0, 10.0]);
        assert_eq!(g_list, vec![1.0, 1.0]);
    }

    #[test]
    fn test_clean_off_artifacts() {
        let mut item = "A=2BM)".to_string();
        let cleaned_item = clean_off_artifacts(&mut item);
        assert_eq!(cleaned_item, "A=2B");

        let mut item = "B->A + 3C_DUP".to_string();
        let cleaned_item = clean_off_artifacts(&mut item);
        assert_eq!(cleaned_item, "B->A + 3C");
    }

    #[test]
    fn test_search_substances() {
        let mut ReactionAnalyzer_instance = ReactionAnalyzer::new();
        let reactions_: Vec<&str> = vec!["A=2BM)", "B->A + 3C_DUP"];
        let reaction = reactions_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.reactions = reaction;
        ReactionAnalyzer_instance.search_substances();
        assert_eq!(ReactionAnalyzer_instance.substances, vec!["A".to_string(), "B".to_string(), "C".to_string()]);
    }
}