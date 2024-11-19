import numpy as np
from skimage.feature import peak_local_max
from scipy.special import gamma
from scipy.optimize import minimize
import re

def GetCarpetSignal(hist : np.ndarray) -> np.ndarray: # Достаёт из полной гистограммы отклика центральный "Ковёр"
    return hist[hist.shape[0] // 2 - 10:hist.shape[0] // 2 + 10, hist.shape[0] // 2 - 10:hist.shape[0] // 2 + 10]

def VEM(hist : np.ndarray) -> np.ndarray: # Конвертирует сигнал в р.ч. в показания LC-преобразователей, затем - обратно в р.ч. Имитируем дискретность реальных данных.
    filt = (hist < 0.5)
    hist[filt] = 0.001
    hist = (1 + np.log(hist/8)/np.log(1.12)) // 1
    hist[hist < 0] = 0
    rel_p = 8*1.12**(hist - 1)
    rel_p[rel_p < 8] = 0

    return rel_p

def VEM2rp(hist : np.ndarray) -> np.ndarray: # Конвертирует сигнал LC-преобразователей в р.ч. 
    
    rel_p = 8*1.12**(hist - 1)
    rel_p[rel_p < 8] = 0

    return rel_p

#Функция ниже определяет положение оси ливня по локальному максимуму сигнала в центральном "Ковре"

def neighbors(matrix, rowNumber, colNumber):
    
    result = []
    for rowAdd in range(-1, 2):
        newRow = rowNumber + rowAdd
        if newRow >= 0 and newRow <= len(matrix)-1:
            for colAdd in range(-1, 2):
                newCol = colNumber + colAdd
                if newCol >= 0 and newCol <= len(matrix)-1:
                    if newCol == colNumber and newRow == rowNumber:
                        continue
                    result.append(matrix[newCol][newRow])
    return np.array(result)

def find_smoothest_max(signal, coords):

    if (len(coords) == 0):
        return np.array([0, 0])

    decrease = np.zeros(len(coords))
    
    for i, coordinate in enumerate(coords):

        amplitude = signal[coordinate[0], coordinate[1]]

        neighbours = neighbors(signal, coordinate[1], coordinate[0])

        decrease[i] = np.sum(neighbours) / (len(neighbours))

    smoothest_max_coord = coords[np.argmax(decrease)]

    return smoothest_max_coord

def get_xy_lm(signal):

    x_cover = np.round(np.arange(-6.65, 6.66, 0.7), 2) # Координаты детекторов Ковра (м)
    y_cover = np.round(np.arange(-6.65, 6.66, 0.7), 2) # Координаты детекторов Ковра (м)

    coordinates = peak_local_max(signal, exclude_border=False, threshold_rel=0.1)
    i, j = find_smoothest_max(signal, coordinates)

    x = x_cover[j]
    y = y_cover[19 - i]

    return [x, y]
    
# Функция ниже методом МНК строит линейную зависимость индексов детекторов с максимальным энерговыделением в каждой строке (столбце) от номера строки (столбца). 
# В роли весов выступают суммарные энерговыделения в соответсувующей строке (столбце).
# На выходе получаем коэффеициенты соответсвующих полиномов kx + b

def fit_coeffs(row_ind : np.ndarray, row_centres : np.ndarray, row_weights : np.ndarray, 
               col_ind : np.ndarray, col_centres : np.ndarray, col_weights : np.ndarray) -> np.ndarray:

    row_polyfit : np.ndarray = np.polynomial.Polynomial.fit(row_ind, row_centres, deg=1, w=row_weights).convert().coef # Коэффициенты аппроксимаций, полученные методом МНК
    col_polyfit : np.ndarray = np.polynomial.Polynomial.fit(col_ind, col_centres, deg=1, w=col_weights).convert().coef # Коэффициенты аппроксимаций, полученные методом МНК

    b_row : float = row_polyfit[0] #Явно объявляем коэффициенты
    k_row : float = row_polyfit[-1] #Явно объявляем коэффициенты

    b_col : float = col_polyfit[0] #Явно объявляем коэффициенты
    k_col : float = col_polyfit[-1] #Явно объявляем коэффициенты

    row_errors : np.ndarray = (row_centres - k_row*row_ind - b_row)**2 #Вычисляем квадраты ошибок аппроксимаций
    col_errors : np.ndarray = (col_centres - k_col*col_ind - b_col)**2 #Вычисляем квадраты ошибок аппроксимаций

    # В случае наличия ошибки на более, чем 5 рядов/столбцов выкидываем худшую точку из аппроксимации (защита от вбросов)
    # Далее проводим аппроксимацию заново

    if np.any(row_errors > 5**2): 
        
        worst_point : int = np.argmax(row_errors)
        
        row_ind : np.ndarray = np.delete(row_ind, worst_point)
        row_centres : np.ndarray = np.delete(row_centres, worst_point)
        row_weights : np.ndarray = np.delete(row_weights, worst_point)

        row_polyfit : np.ndarray = np.polynomial.Polynomial.fit(row_ind, row_centres, 1, w=row_weights).convert().coef 

        b_row : float = row_polyfit[0]
        k_row : float = row_polyfit[-1]

    if np.any(col_errors > 5**2):
        
        worst_point : int = np.argmax(col_errors)
        
        col_ind : np.ndarray = np.delete(col_ind, worst_point)
        col_centres : np.ndarray = np.delete(col_centres, worst_point)
        col_weights : np.ndarray = np.delete(col_weights, worst_point)

        col_polyfit = np.polynomial.Polynomial.fit(col_ind, col_centres, 1, w=col_weights).convert().coef

        b_col : float = col_polyfit[0]
        k_col : float = col_polyfit[-1]

    return np.array([[b_row, k_row], [b_col, k_col]])

# Функция ниже определяет положение оси ливня (номер детектора) по отклику Ковра, подаваемому на вход

def get_xy_fits(hist : np.ndarray) -> list[float]:

    # Первым делом вычисляется количество ненулевых строк и столбцов. Методу необходимо минимум две непустые строки и столбца.
    # В противном случае в качестве оси указывается центр Ковра

    non_zero_cols : int = np.count_nonzero(np.count_nonzero(hist, axis=0))
    non_zero_rows : int = np.count_nonzero(np.count_nonzero(hist, axis=1))
    
    if ((non_zero_cols >= 2) and (non_zero_rows >= 2)):
    
        n_rows : int = hist.shape[0]
        n_cols : int = hist.shape[1]
    
        row_ind : np.ndarray = np.arange(0, n_rows) # Номера строк
        col_ind : np.ndarray = np.arange(0, n_cols) # Номера столбцов
    
        row_weights : np.ndarray = np.zeros(n_rows)
        row_centres : np.ndarray = np.zeros(n_rows)
        
        col_weights : np.ndarray = np.zeros(n_cols)
        col_centres : np.ndarray = np.zeros(n_cols)

        # Для каждой строки (столбца) определяется индекс k детектора с максимальным энерговыделением. 
        # Таким образом получаем точку на плоскости: х - номер строки (столбца), y - k
        # Вес точки - суммарное энерговыделение в строке (столбце)
    
        for i in range(n_rows):
            
            row : np.ndarray = hist[i]
            
            weight : float = np.sum(row)

            maxima : int = np.argwhere(row == row.max())
            centre : float = np.mean(maxima)
            
            row_weights[i] = weight
            row_centres[i] = centre
            
        for i in range(n_cols):
    
            col : np.ndarray = hist[:, i]
            
            weight : float = np.sum(col)
            
            maxima : int = np.argwhere(col == col.max())
            centre : float = np.mean(maxima)
            
            col_weights[i] = weight
            col_centres[i] = centre
    
        coeffs : np.ndarray = fit_coeffs(row_ind, row_centres, row_weights, col_ind, col_centres, col_weights) # Коэффициенты полученных аппроксимаций

        b_row : float = coeffs[0, 0] # Явно объявляем коэффициенты
        k_row : float = coeffs[0, 1] # Явно объявляем коэффициенты

        b_col : float = coeffs[1, 0] # Явно объявляем коэффициенты
        k_col : float = coeffs[1, 1] # Явно объявляем коэффициенты

        if (k_col*k_row == 1): # Случай k_col*k_row == 1 означает, что линии совпадают. В таком случае осью ШАЛ считается середина полученной линии
            col0 : float = np.mean(col_ind[np.nonzero(col_weights)])
            row0 : float = k_col*col0 + b_col
        else: # Поиск координат пересечения двух прямых на плоскости
            col0 : float = (k_row*b_col + b_row) / (1-k_row*k_col)
            row0 : float = k_col*col0 + b_col
        
        return [row0, col0]
    
    else:
        
        return [9.5, 9.5]

# Функция ниже определяет положение оси ливня (в метрах) по отклику Ковра, подаваемому на вход

def get_xy(signal : np.ndarray) -> list[float]: 

    i, j = get_xy_fits(signal) # Определение положения оси по номерам детекторов

    x : float = 0.7*j - 6.65 # Положение оси ливня в метрах
    y : float = 6.65 - 0.7*i # Положение оси ливня в метрах

    return [x, y]

# Функция ниже аналитически определяет зенитный угол прихода ШАЛ в приближении плоского фронта

def get_PFA_theta(t : np.ndarray) -> tuple:

    x : np.ndarray = np.array([-21.29, -21.63, 17.5, 21.16]) # Координаты выносных пунктов Ковра (м)
    y : np.ndarray = np.array([21.68, -21.13, 9.34, -21.05]) # Координаты выносных пунктов Ковра (м)
    z : float = 1700 # Примерная высота над уровнем моря для Ковра (м)
    
    c_norm : float = 0.3 # Скорость света, нормированная на 10**(-9), т.к. позже будем домножать на наносекунды
    
    x_sq : np.ndarray = x**2 # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y_sq : np.ndarray = y**2 # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    xy : np.ndarray = x*y # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    
    x_mean : float = x.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y_mean : float = y.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    x_mean_sq : float = x_mean**2 # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y_mean_sq : float = y_mean**2 # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    x_sq_mean : float = x_sq.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y_sq_mean : float = y_sq.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    xy_mean : float = xy.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    xy_mean_sq : float = xy_mean**2 # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    
    t_mean : float = t.mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    xt_mean : float = (x*t).mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    yt_mean : float = (y*t).mean() # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)

    xt_dif : float = xt_mean - x_mean*t_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    x2t_dif : float = x_sq_mean*t_mean - x_mean*xt_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    x2x_dif : float = x_mean_sq - x_sq_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)

    yt_dif : float = yt_mean - y_mean*t_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y2t_dif : float = y_sq_mean*t_mean - y_mean*yt_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)
    y2y_dif : float = y_mean_sq - y_sq_mean # Данные необходимые для аналитического расчёта направления прихода ШАЛ (в приближении плоского фронта)

    nx : float = (xy_mean*yt_dif+x_mean*y2t_dif+xt_mean*y2y_dif)/(x_sq_mean*y_mean_sq+x_mean_sq*y_sq_mean-2*x_mean*y_mean*xy_mean+xy_mean_sq-x_sq_mean*y_sq_mean)*c_norm
    ny : float = (xy_mean*xt_dif+y_mean*x2t_dif+yt_mean*x2x_dif)/(x_sq_mean*y_mean_sq+x_mean_sq*y_sq_mean-2*x_mean*y_mean*xy_mean+xy_mean_sq-x_sq_mean*y_sq_mean)*c_norm

    nz : float = np.sqrt(1-nx**2-ny**2) # Иногда 1-nx**2-ny**2 оказывается меньше 0, что означает что в данном случае приближение плоского фронта не работает 

    theta : float = np.arccos(nz)
    
    if np.isnan(theta):
        return (-1, -1, -1, -1)
    else:
        return (np.degrees(theta), nx, ny, nz)

# Функция ниже вычисляет хи-квадрат (отклонение данных от модифицированной модели плоского фронта)

def Chi_sq(params : tuple, x0 : float, y0 : float, t: np.ndarray) -> float:

    x : np.ndarray = np.array([-21.29, -21.63, 17.5, 21.16]) # Координаты выносных пунктов Ковра (м)
    y : np.ndarray = np.array([21.68, -21.13, 9.34, -21.05]) # Координаты выносных пунктов Ковра (м)
    z : float = 1700

    sigma_to : float = 2.6 # Параметры модифицированной модели плоского фронта ШАЛ, учитывающей ослабление ливня по мере удаления от оси
    b : float = 1.5 # Параметры модифицированной модели плоского фронта ШАЛ, учитывающей ослабление ливня по мере удаления от оси
    r_t : int = 30 # Параметры модифицированной модели плоского фронта ШАЛ, учитывающей ослабление ливня по мере удаления от оси

    c_norm : float = 0.3 
    
    nx, ny, nz, t0 = params

    chi_sq : float = 0

    n : int = len(t)
    
    for i in range(n):
        r_i : float = np.sqrt((x[i]-x0)**2+(y[i]-y0)**2)
        sigma_i : float = sigma_to*(1+r_i/r_t)**b
        w_i : float = 1/(sigma_i**2)

        chi_sq += w_i*(nx*x[i]+ny*y[i]*nz*z-c_norm*(t[i]-t0))**2
        
    return chi_sq

def constr(pars : tuple) -> float: # Ограничение на сумму nx, ny, nz. Оно используется при минимизации хи-квадрата
    nx, ny, nz, t0 = pars

    return nx**2+ny**2+nz**2-1

def get_PFAWTC_theta(t : np.ndarray, x0 : float, y0: float) -> np.ndarray: # Данная функция выдаёт тета, фи ШАЛ в модифицированном приближении плоского фронта

    constraint = {'type': 'eq', 'fun': constr} # Задание ограничения выше в формате scipy

    PFA_result : tuple = get_PFA_theta(t)

    if (PFA_result == (-1, -1, -1, -1)):
        return [-1, -1]

    initial_guess : np.ndarray = np.append(PFA_result[1:], t.min())

    result = minimize(Chi_sq, initial_guess, args=(float(x0), float(y0), t), constraints=constraint, tol=1e-6)

    theta : float = np.arccos(result.x[2]) #Расчёт тета по параметру nz
    phi : float = (1 - np.sign(result.x[0]))*np.pi/2 + (1 + np.sign(result.x[0]))*(1 - np.sign(result.x[1]))*np.pi/2 + np.arctan(result.x[1]/result.x[0]) # Расчёт фи с учётом четверти
    phi_moved : float = np.radians(242) - phi # Здесь учитывается то, что Ковёр повернут относительно вектора, указывающего на север
    if phi_moved < 0:
        phi_moved += 2*np.pi

    return np.array([theta, phi_moved])

# Данная функция рассчитывает поверхностную плотность частиц в детекторах Ковра с учётом переходного эффекта
def get_rho(signal : np.ndarray, time_array : np.ndarray = None, angles : np.ndarray = None) -> list[np.ndarray]:

    x_cover = np.round(np.arange(-6.65, 6.66, 0.7), 2) # Координаты детекторов Ковра (м)
    y_cover = np.round(np.arange(-6.65, 6.66, 0.7), 2) # Координаты детекторов Ковра (м)

    x0, y0 = get_xy(signal) # Сначала вычисляется положение оси ливня
    
    if (np.any(angles)):
    	phi, theta = angles
    else:
    	theta, phi = get_PFAWTC_theta(time_array, x0, y0) # С учётом вычисленного положения оси получаем тета, фи

    rho : list = []
    r : list = [] 

    for i in range(len(y_cover)):
        for j in range(len(x_cover)):

            x : float = x_cover[j]
            y : float = y_cover[19 - i]
            
            # Вычисление расстояния до детектора в плоскости ливня
            r_ij : float = np.sqrt(((x-x0)*np.sin(theta)*np.sin(phi) - (y-y0)*np.cos(theta)*np.cos(phi))**2 + ((x-x0)*np.cos(theta))**2 + ((y-y0)*np.cos(theta))**2)
            
            k_ij : float = 1 + 7.503 / (1.636 + r_ij**1.474) # Переходный эффект крыши здания, в котором находится Ковёр
            rho_ij : float = signal[i, j] / (0.49 * k_ij)

            if (r_ij > 0) and (rho_ij > 0):

                rho.append(rho_ij)
                r.append(r_ij)

    r : np.ndarray = np.array(r)
    rho : np.ndarray = np.array(rho)
    ind : np.ndarray = np.argsort(r)

    return [r[ind], rho[ind]]

def NKG(x : float, s: float, Ne : int) -> float:
    return (Ne / 95**2) * (gamma(4.5 - s)/(2*np.pi*gamma(s)*gamma(4.5 - 2*s))) * (x / 95)**(s - 2) * (1 + (x/95))**(s - 4.5)

def NKG_Vik(x : float, s: float, Ne : int) -> float: 
    return ((Ne/(95**2))*((5.803**(-(s-1.26)**2))/2.26)*(x/95)**(s-2)*(1+x/95)**(s-4.5))

def num_sort(test_string : str) -> list:
    return list(map(int, re.findall(r'\d+', test_string)))[0]
    

    
# Функция "мексиканская шляпа" также может использоваться для определения положения оси ШАЛ
# Истинному положению соответсвует максимуму этой функции
# Для простоты hat_func ниже возвращает значение функции, домноженное на -1, чтобы свести задачу поиска максимума к поиску минимума

def hat_func(params : list[float], hist : np.ndarray) -> float:

    y, x = params
    
    n_rows : int = hist.shape[0]
    n_cols : int = hist.shape[1]

    f : float = 0

    for i in range(n_rows):
        for j in range(n_cols):
            
            row : int = i
            column : int = j

            t_ij_sq : float = (row - y)**2 + (column - x)**2

            entry : float = hist[row, column]*(2 - t_ij_sq)*np.exp(-t_ij_sq/2)

            f += entry

    return -f
    
# Эта функция осуществляет поиск минимума функции "мексиканской шляпы" методом имитации отжига
# Метод стохастический, т.е. случайный. На каждом шаге случайно генерируется х, у
# Каждой точке сопоставляется значение "энергии" (в нашем случае функции "шляпы")
# Алгоритм находит точку минимуму "энергии"

def get_annealing_xy(hist : np.ndarray) -> list[float]:

    T_max : int = 25 # Диапазон "температур" для имитации отжига
    T_min : float = 0.01 # Диапазон "температур" для имитации отжига

    n_rows : int = hist.shape[0]
    n_cols : int = hist.shape[1]

    T : float = T_max
    y : float = 9.5
    x : float = 9.5
    E : float = hat_func([y, x], hist)

    for i in range(1, 1_000):
        
        x_rand : float = np.random.uniform(-0.5, n_cols - 0.5)
        y_rand : float = np.random.uniform(-0.5, n_rows - 0.5)
        E_rand : float = hat_func([y_rand, x_rand], hist)
        deltaE : float = E_rand - E

        if deltaE < 0:
            
            y = y_rand
            x = x_rand
            E = E_rand
            
        else:
            
            prob : float = np.exp(-deltaE/T)
            change_prob : float = np.random.uniform()

            if change_prob <= prob:

                y = y_rand
                x = x_rand
                E = E_rand

        T = T_max/np.log(1 + i)

        if T <= T_min: break

    return [y, x]  
