import numpy as np
from collections.abc import Callable
from scipy.optimize import OptimizeResult
from abc import ABC, abstractmethod

#Librería para encontrar raíces de funciones de una variable. Creada por Fabrizio Benavidez 12/03/2026
#Library to find roots of one variable functions. By Fabrizio Benavidez
#This .py is a preliminary idea about how to find roots of a function of one variable. Therefore there will be modifications.


#Decoradores para mi clase

# def Bolzano(func):
#     def wrapper(f, xl, xu, *args):
#         if f(xl)*f(xu) > 0:
#             message = f'There is not a change of sign between {xl} and {xu}'
#             raise ValueError(message)
#         return func(f, xl, xu, *args)
#     return wrapper

#Debe ser más conveniente este
def Bolzano(method): #AL ser un método es mejor usar method como nombre. Ahora este decorador es universal
    #Al usar self ya puedo entrar a todas las instancias de mi clase
    def wrapper(self, *args):
        # Aquí validas usando self.f, self.xl, etc.
        if np.sign(self.f(self.xl)) == np.sign(self.f(self.xu)):
            message = f'There is not a change of sign between {self.xl} and {self.xu}'
            raise ValueError(message)
        return method(self, *args)
    return wrapper

def ByZero(method):
    def wrapper(self, *args):
        if np.allclose(self._df(self.xi), 0, atol = 1e-12):
            message = f'The derivative in {self.xi} tends to 0'
            raise ZeroDivisionError(message)
        return method(self, *args)
    return wrapper
        



#Clase abstracta
class MetodoNumericoCerrado(ABC):

    @abstractmethod
    def __init__(self, f : Callable, xl, xu):
        self.f = f
        self.xl = xl
        self.xu = xu
        
    
    @abstractmethod
    def resolver(self):...

class MetodoNumericoAbierto(ABC):

    @abstractmethod
    def __init__(self, f : Callable, xi):
        self.f = f
        self.xi = xi
        
    
    @abstractmethod
    def resolver(self):...


#Mi método de raíces
class Biseccion(MetodoNumericoCerrado):
    es = None
    max_iter = 100

    def __init__(self, f, xl, xu):
        super().__init__(f, xl, xu)

    @Bolzano
    def resolver(self):
        if not callable(self.f):
            message = 'f must be a Callable'
            raise ValueError(message)
    

        if self.es is None:
            self.es = 1e-8
        else:
            self.es = self.es

        success = True
        iter = 0
        xr = self.xu
        ea = 100

        while self.es <= ea:
            if iter >= self.max_iter:
                cause = RuntimeWarning('max iter reached')
                success = False
                return OptimizeResult(xr = xr, success = success, i = iter, ea = ea, message = cause)
            xv = xr
            xr = (self.xl + self.xu)/2
            if xr != 0:
                ea = np.abs((xr - xv)/xr)*100
            iter += 1
            if np.allclose(self.f(xr), 0, atol = 1e-12):
                break
            proof = self.f(self.xl)*self.f(xr)
            if proof < 0:
                self.xu = xr
            else:
                self.xl = xr
        
        return OptimizeResult(xr = xr,
                        success = success,
                        i = iter,
                        ea = ea)


class NewtonRaphson(MetodoNumericoAbierto):
    es = None
    max_iter = 100

    def __init__(self, f, xi, df : Callable) -> None:
        super().__init__(f, xi)
        self._df = df

    @ByZero
    def resolver(self):
        if not (callable(self.f) and callable(self._df)):
            message = 'f and df must be callables'
            raise ValueError(message)

        if self.es is None:
            self.es = 1e-8
        else:
            self.es = self.es

        success = True
        iter = 0
        ea = 100 

        while self.es <= ea:
            if iter >= self.max_iter:
                cause = RuntimeWarning('max iter reached')
                success = False
                return OptimizeResult(xr = xr, success = success, i = iter, ea = ea, message = cause)
            xr = self.xi - self.f(self.xi)/self._df(self.xi)
            if xr != 0:
                ea = np.abs((xr - self.xi)/xr)*100
            iter += 1
            self.xi = xr

        return OptimizeResult(xr = xr,
                        success = success,
                        i = iter,
                        ea = ea)
        





if __name__ == '__main__':

    def func(x):
        return -np.pi + np.exp(-0.14586*x) 
    
    def df_func(x):
        return -0.14586*np.exp(-0.14586*x)

    def func2(x):
        return (x-4)**2 - 2
    
    def df_func2(x):
        return 2*(x-4)


    xlower = -10
    x_initial = 2
    xupper = 0

    root_func = Biseccion(func, xlower, xupper)
    root_func_NR = NewtonRaphson(func2, x_initial, df=df_func2)
    print('Método de bisección: \n', root_func.resolver())
    print('Método de Newton-Raphson: \n', root_func_NR.resolver())
