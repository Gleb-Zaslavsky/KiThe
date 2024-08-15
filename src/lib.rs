/// ru 
/// Модуль берет на вход  вектор уравнений реации, заданных в виде String и выдает следующие данные: 
/// 1) стехиометрическую матрицу заданную в виде вектора векторов 
/// 2) вектор веществ 
/// 3) вектор векторов стехиометрических коэффициентов реагентов в каждой реакции 
/// 4) то же для продуктов
/// В процессе производится избавление от артефактов парсинга уравнений из баз данных
/// Примечание: 
/// 1) уравнения реакции содержат последние символы '_dup' или '_DUP' (то есть дублирующие) - эти 
/// символы удаляются как артефакты парсинга
/// 2) для работы с эмпирическими реакциями у которых 
///  код возвращает следующие структуры данных: матрица стехеометрических коэффициентов,
///  матрица коэффициентов прямых реакций и матрица коэффициентов обратных реакций, матрица степеней концентраций для кинетической функции, G_matrix,
/// как правило степени концентраций в кинетической функции совпадают со стехеометрическими коэффициентами веществ в реакциии, однако,
/// для эмпирических реакций они могут и отличаться от стехеометрических коэффициентов.
///  Предусмотрена вохможность прямого указания этих коэффициентов в уравнениях реакции в виде "степени"
/// после формулы вещества (например A**0.3 ) такие степени записываются в  G_matrix вместо стехеометрических коэффициентов
/// в противном случае записывается стехеометрический коэффициент
/// ----------------------------------------------------------------
/// eng
/// The module takes as input a vector of reaction equations specified as a vector of String and produces the following data:
 /// 1) a stoichiometric matrix specified as a vector of vectors 
 /// 2) a vector of substances 
 /// 3) a vector of vectors of stoichiometric coefficients of reactants in each reaction 
 /// 4) the same for products
/// The process involves getting rid of artifacts of parsing equations from databases
/// 
/// Note: 
/// 1) reaction equations contain the last characters '_dup' or '_DUP' (that is, duplicate ones) - these 
/// characters are removed as parsing artifacts
/// 2) for working with empirical reactions in which  the code returns the following data structures: matrix of stoicheometric coefficients,
/// matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the 
/// kinetic function, G_matrix. As a rule, the degrees of concentration in the kinetic function coincide with the stoicheometric coefficients of 
/// the substances in the reaction; however, for empirical reactions they may differ from the stoicheometric coefficients.
/// It is possible to directly indicate these coefficients in the reaction equations in the form of a “degree” after the formula of a substance 
/// (for example A**0.3), such degrees are written in G_matrix instead of stoicheometric coefficients otherwise, the stoicheometric coefficient is written
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

pub mod reaction_analyzer;
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
/// processing of user-chosen reactions.
///  So you define what reactions you need by using the constructor of mechanism from the mechfinder_api module or
/// manually with help of kinetics_lib_api module. Now you ca
pub mod User_reactions;
// Basis functionality to search in reaction library
pub mod kinetics_lib_api;
/// ru
///  в проекте используются шорткаты для названий реакций в стиле
/// "сокращенное название библиотеки реакции"+"номер реакции в библиотеке",
///  например C_10 означает реакция '10' в библиотеке CANTERA. 
/// То же означает C10 и Сantera_10
/// eng
/// The project uses shortcuts for reaction names in the style 
/// of "short name of the reaction library" + "reaction number in the library" 
/// for example C_10 means reaction '10' in the CANTERA library. 
/// The same means C10 and Сantera_10
pub mod parsetask;