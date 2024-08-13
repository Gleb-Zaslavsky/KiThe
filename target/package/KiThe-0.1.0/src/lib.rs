
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

pub mod react_vec_parse;
 /// ru 
 /// Модуль снабжен библиотекой кинетических параметров химических реакций, полученной в результате парсинга общедоступныхбаз данных 
 /// Модуль берет на вход название библиотеки и вектор веществ а затем выдает следующие данные:
 /// 1) все реакции исходных веществ между собой, и всех их возможных продуктов между собой. 
 /// 2) HashMap с кинетическими данными всех найденных реакций
/// ----------------------------------------------------------------
 /// eng
 /// The module is equipped with a library of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases 
 /// The module takes as input the name of the library and the vector of substances and then produces the following data:
 /// 1) all reactions of starting substances with each other, and all their possible products with each other. 
 /// 2) HashMap with kinetic data of all found reactions
pub mod mechfinder_api;

