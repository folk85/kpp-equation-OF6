2021-02-06

- `myBturbFOAM`: Add correlation curves in `Rx` vector. Calculate integral length scale. 

2021-02-05

- `myBturbFOAM`: Integrate FFT for velocity to get spectral kinetic energy distribution `Ek`. Use three types of EWM with alpha 1e-3, 1e-5 and `dt/100`. Add calculation of autocorrelation function `Rx`.

2021-01-14

- Set it 

2020-12-24

- Add Stochastic Burgers equation with  initial spectrum of velocity by Ornstein-Uhlenbeck.

2020-07-08

- Optimize the Random field generation. Replace dimensional vars to scalars. It makes good performance for time depended Random fields 
- Add Random telegraph Noise 

2020-04-29

- Set Normalized equation with Dynamic Orstein-Uhlenbbeck named (sdeScheme == "OrsteinTimeN"). 
- Define  mean velocity at every timeStep ($<V> = \int_{0}^{L} U*(1-U) dx$) with Dynamic Orstein-Uhlenbbeck named (sdeScheme == "OrsteinTimeM"). 

2020-04-22

- Add Dynamic White Noise field.

2020-04-18

- Fix dynamic OU process according the new recommedations. Dispersion of the system equals  SQRT(\tau * \theta).

- Refuse the velocity shifting to set Mean Value to Zero. This operation was used for stochastic field by space only. 

- Add the IO infomation of theta, tau, barVel and SDE_Scheme

2020-04-14

Set dynamical stochastic fields, which changing during the time. 

- OrsteinTime -- generate 2D Orstein-Uhlenbeck stochastic field by scheme [Krasheninnikov]. It considers exponentional delay by time  and space
- WhiteNoiseTime -- generate new white noise field at every timestep.

2020-04-10

Add option `sdeScheme` to choose the Stochastic field. Now there are two options:

- Ornstein  -- generate Ornstein-Uhlenbeck stochastic field
- WhiteNoise -- generate White noise velocity field

Currently both gen method work in one-D X-space. 

By Default set Stochastic field Orstein-Uhlenbeck along the X-space as an initial conditions
It works with the old version of cases. If the variable "sdeScheme" isn't set in transportProperties,
then the program uses the O-U process. 



В уравнение со сверткой myKppOUNeFoam Добавили код с выделением среднего значения винеровского процесса в области пламени. Проверяли возможность применения такого преобразования для упрощения решения при больших амплитудах. Сейчас код закомментирован. 

2020-03-16

Добавлил систему для экспоненциального белого шума (Процесс Орнштейна-Уленбека). В модели найдена ошибка, которую раньше не замечали . случайную величину ненужно умножать на sqrt(dx). Это никак не проверить на уже сгенерированном поле. Не влияет на корреляционную функцию. Значение можно вынести в значение амплитуды \delta.

2020-03-13

- Разнесли приложения по директориям. назвали разными  названиями. Сейчас два варианта: расчет с белым шумом  и расчет с белым шумом в альтернативной записи
- Исправлена ошибка при вычислении экспоненциальной переменной. После решения уравнения проводится коррекция решения .


2020-03-12

Введена альтернативная запись через свертку оператора лапласса и конвективного члена. В новой постановке задача становится более гладкой. Проверим ка это влияет на результаты расчетов с большими амплитудами возмущений.

Добавлена переменная для анализа div(U). Источники  записаны линеаризованном виде. Белый шум записан в виде N(0,1)/Sqrt(h)

- 2020-02-19: 
    
    - Добавлена переменная для стохастического поля. Скорость вычисляется как произведение velInit на амплитуду.
    - Скорректировано уравнение по знакам, чтобы синхронизироваться с Посвянским. Теперь 
    $${  }$$



Определим уравнение КПП со стохастическим параметром