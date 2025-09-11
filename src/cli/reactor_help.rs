pub const REACTOR_RU_HELPER:&'static str = 

"
                                Общие соображения. \n
При использовании программы следует помнить, что все величины необходимо указывать в системе СИ. Это стоит иметь в виду и \n
тщательно проверять размерность вводимых физических величин! \n
1) Инструкции файла задания состоят из заголовков – строк, содержащих одно слово, под которым размещаются пары ключ - значения,\n
всегда разделенные знаком двоеточия. Заголовок всегда означает некоторую логически обособленную группу параметров, ключи, \n
именуемые полями – имена параметров, а значения – это собственно значения параметров. Параметры могут быть целочисленными (100),\n
числами с плавающей точкой (13.054, 1e-4), строковыми, векторными ([1.0, 2.0, 3.0]), булевыми, а также опциональными, это \n
означает, что указание вместо параметра None это легитимная возможность – в этом случае будет задействовано дефолтное \n
значение. Такие опциональные параметры оборачиваются в Some(..) например Some(true), Some(1e-3) и т.д. \n
2) Заголовки и поля могут быть опциональными или обязательным. Неуказание обязательного заголовка или поля (естественно с \n 
некоторым значением) вызовет вылет программы.  Неуказание опционального заголовка или поля означает, что либо «под капотом» \n 
будет выбрано значение по умолчанию, либо некоторые дополнительные опции не будут задействованы. \n

Итак, опишем ниже все заголовки, поля и их значения с указанием обязательности и типичных значений.\n

                                Задание условий процесса \n			
process_conditions – обязательный заголовок. Все его поля, тоже являются обязательными если иное не указано прямо.\n
problem_name: Some(%string%), \n
problem_description: Some(%string%) – поля содержащие название и описание задачи. Тип значения опциональная строка, что \n
означает, что это значение может быть None. Ни на что не влияет. \n
substances: %substance_1%, %substance_2%,… - поле с перечнем веществ. \n
        t0: 0.0 – начальное значение аргумента; \n
        t_end: 1.0 – конечное значение аргумента; \n
Поскольку характерная длина задачи указана ниже, а в программе используется масштабированные величины, реальный диапазон \n
 значений аргумента будет от L*t0 до L*t_end поэтому логично будет просто принять их как 0 и 1.\n
        n_steps: %integer% поле со значением типа целое для указания количества  \n
шагов дискретизации (см. ниже);\n
        arg:%string% – необязательное поле со значением строкового типа, содержит название аргумента, если поле со  \n
        значением не задано, аргумент по умолчанию будет назван “x”. \n
        Tm: %float% - значение температуры при которой вычисляются значения термодинамических величин – среднее значение \n
         температуры во фронте горения, К.\n
        L: %float% - оценка характерного размера задачи, м\n
В программе применяется следующая формула для перехода к безразмерной температуре Teta = (T – dT)/T_scale, обычно dT \n
полагают равной или несколько меньше начальному условию для температуры T0, а T_scale полагают равным приросту температуры \n
в системе, то есть разнице между T-max и T0  \n
        dT: %float%  , К \n
        T_scale: %float%  ,К \n
        P: %float% - давление в системе, Па \n
        Cp: %float% - теплоемкость, Дж/кг*К \n
        Lambda: %float% - коэффициент теплопроводности, Вт/м*К \n
        m: %float% - массовая скорость горения, кг/м2*с \n
        M: %float% - необязательное поле, если его не указать будет вычислено автоматически, кг/моль \n
                                n_steps & L \n
Под характерной длиной для данной задачи будем далее иметь в виду такое пространственное расстояние, на котором все \n
химические реакции заведомо закончились, а градиенты температуры и концентраций сделались близкими к нулю. Если L выбрать \n
много большим, чем характерное расстояние, то при относительно небольшом количестве точек (n_steps) на фронт химических \n
реакций будет приходиться малое количество точек сетки, солвер не сможет разрешить больших градиентов и задача вылетит с \n
ошибкой. С другой стороны, задание слишком большого количества точек будет затратно в вычислительном смысле.  \n 
Если же характерное расстояние выбрать слишком малым, то решение может отсутствовать, действительно, условия на правой  \n
границе требуют обращения в нуль градиентов концентраций и температуры. Таким образом, L следует выбрать несколько большим \n
чем характерная длина задачи, тогда это не потребует чрезмерно большого количества точек. Предположим, мы выбрали L \n
разумно, как же мы должны выбрать n-steps? Если решено отказаться от использования адаптивных сеток (adaptive: None), \n
то достаточно 100-200 точек. В случае если решено задействовать адаптивные сетки, можно начать с 10-20 точек, в этом \n
случае алгоритм перерасчета сеток сам добавит достаточно точек в тех областях сетки, где наблюдаются большие градиенты \n
параметров. По поводу выбора алгоритмов grid refinement и параметров этих алгоритмов см. ниже

                                Постпроцессинг.\n
Программа предлагает несколько опций что делать с решением. Для того, чтобы задействовать опции постпроцессинга необходимо \n
объявить в файле задания опциональный заголовок postprocessing, под которым уже разместить желаемые инструкции построцессинга. \n
Итак, ниже вы можете указать следующие поля (все они без исключения опциональны) со значениями:\n
* gnuplot:true вызывает gnuplot, бесплатную программу для построения двух- и трёхмерных графиков, для построения решения. \n
 Очевидно, что эта программа должна быть установлена на компьютере пользователя, иначе gnuplot:true приведёт к ошибке. Как \n
 уже упоминалось, это бесплатная программа, которую можно скачать из Интернета и установить. Также необходимо убедиться, что \n
она добавлена в Path. \n
*plot:true вызывает нативную библиотеку Rust для построения графиков. Не требует внешних зависимостей. \n
* save:true команда сохранения решения в файл txt \n
* solve_to_csv:true  команда сохранения решения в файл csv формата. \n
Если не задано поле filename со значением %имя_которое_вы_выбрали% файлы txt и csv будут названы result. Поэтому опциональное \n
 поле filename отвечает за имена создаваемых файлов. \n
* return_to_dimension – опциональное поле с булевым значением - означает переходить ли от безразмерных переменных к размерным.\n
 Значение по умолчанию true – то есть если поле вообще не задать будет установлено именно это значение. Однако, иногда \n
 полезно посмотреть на безразмерные переменные в этом случае можно указать return_to_dimension: false \n



";
pub const REACTOR_ENG_HELPER:&'static str =
"
                                General considerations. \n
When using the program, remember that all values must be specified in the SI system. Keep this in mind and carefully check the \n
dimensions of the physical quantities you enter!\n
1) The task file instructions consist of headers—lines containing a single word, followed by key-value pairs, always \n
separated by a colon. A header always denotes a logically separate group of parameters, keys called fields – parameter names, \n
and values – the actual parameter values. Parameters can be integers (100), floating point numbers (13.054, 1e-4), strings, \n
vectors ([1.0, 2.0, 3.0]), Booleans, and optional, which means that specifying None instead of a parameter is a legitimate \n
option—in this case, the default value will be used. Such optional parameters are wrapped in Some(..), for example, \n
Some(true), Some(1e-3), etc.
2) Headers and fields can be optional or mandatory. Failure to specify a mandatory header or field (naturally with some value) \n
will cause the program to crash.  Failure to specify an optional header or field means that either the default value will be \n
selected “under the hood” or some additional options will not be used. \n
So, let's describe all the headers, fields, and their values below, indicating their mandatory nature and typical values.\n

                                Setting process conditions\n
process_conditions – mandatory header. All its fields are also mandatory unless otherwise specified.\n
problem_name: Some(%string%),\n
problem_description: Some(%string%) – fields containing the name and description of the task. The value type is an optional \n
string, which means that this value can be None. It has no effect. \n
substances: %substance_1%, %substance_2%,… – a field with a list of substances.\n
        t0: 0.0 – the initial value of the argument; \n
        t_end: 1.0 – the final value of the argument;\n
Since the characteristic length of the task is specified below and the program uses scaled values, the actual range of \n
argument values will be from L*t0 to L*t_end, so it would be logical to simply accept them as 0 and 1.\n
        n_steps: %integer% field with an integer value to specify the number of discretization steps (see below);\n
        arg:%string% – optional field with a string type value, contains the name of the argument; if the field with the value \n
         is not specified, the default argument will be named “x”. \n
        Tm: %float% - the temperature at which the thermodynamic quantities are calculated – the average temperature at the \n
         combustion front, K. \n
        L: %float% - an estimate of the characteristic size of the problem, m  \n
The program uses the following formula to convert to dimensionless temperature Teta = (T – dT)/T_scale, where dT is usually   \n
assumed to be equal to or slightly less than the initial condition for temperature T0, and T_scale is assumed to be equal to  \n
the temperature increase in the system, i.e., the difference between T-max and T0   \n
        dT: %float%  , K  \n
        T_scale: %float%  ,K  \n
        P: %float% - pressure in the system, Pa \n
        Cp: %float% - heat capacity, J/kg*K \n
        Lambda: %float% - thermal conductivity coefficient, W/m*K \n
        m: %float% - mass combustion rate, kg/m2*s \n
        M: %float% - optional field, if not specified, it will be calculated automatically, kg/mol \n
                                n_steps & L\n
The characteristic length for this problem will be the spatial distance at which all chemical reactions have clearly ended  \n
and the temperature and concentration gradients have become close to zero. If L is chosen to be much larger than the  \n
characteristic distance, then with a relatively small number of points (n_steps), there will be a small number of grid points \n
on the front of the chemical reactions, the solver will not be able to resolve large gradients, and the task will fail with an \n
error. On the other hand, setting too many points will be computationally expensive.  If the characteristic distance is \n
chosen too small, then there may be no solution, since the conditions at the right boundary require the concentration and \n
temperature gradients to be zero. Thus, L should be chosen slightly larger than the characteristic length of the problem,  \n
then it will not require an excessive number of points. Suppose we have chosen L reasonably, how should we choose n-steps? \n
If it is decided to abandon the use of adaptive grids (adaptive: None), then 100-200 points are sufficient. If you decide \n
to use adaptive grids, you can start with 10-20 points, in which case the grid recalculation algorithm will itself add \n
enough points in those areas of the grid where large parameter gradients are observed. For information on selecting grid \n
refinement algorithms and the parameters of these algorithms, see below. \n

                                Postprocessing. \n
The program offers several options for what to do with the solution. To enable postprocessing options, you must declare \n
an optional postprocessing header in the task file, under which you can place the desired postprocessing instructions. \n
Below, you can specify the following fields (all of which are optional) with values: \n
* gnuplot:true calls gnuplot, a free program for constructing two- and three-dimensional graphs, to construct the solution. \n
Obviously, this program must be installed on the user's computer, otherwise gnuplot:true will result in an error. As already \n
mentioned, this is a free program that can be downloaded from the Internet and installed. You also need to make sure that \n
it is added to the Path. \n
*plot:true calls the native Rust library for plotting graphs. No external dependencies are required.  \n
* save:true command saves the solution to a txt file  \n
* solve_to_csv:true  command saves the solution to a csv file. \n
If the filename field is not set to %name_you_chose%, the txt and csv files will be named result. Therefore, the optional \n
filename field is responsible for the names of the files being created. \n

";