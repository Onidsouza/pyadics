from fractions import Fraction

_PADIC_PRECISION = 64  # vamos guardar 10 casas
_MAX_PADIC_EXPONENT = 16  # val maxima é 10, logo a mínima é -9
_MIN_PADIC_EXPONENT = 1 - _MAX_PADIC_EXPONENT
_DISPLAY_CAP = 10  # quando printado, vamos mostra só 10 casas
_USE_UNICODE = True  # imprime ∞ como unicode na tela?

_PRECISION_RANGE = range(0, _PADIC_PRECISION)
_DISPLAY_RANGE = range(0, _DISPLAY_CAP)


def _invmod(a, p, power):  # encontrar um jeito mais rápido usando coeffs
    # usa Hensel's lemma
    c = pow(a, -1, p)
    if power == 1:
        return c
    else:
        b = _invmod(a, p, power-1)
        return (b - (a*b - 1)*c) % p**power


class PAdicFloat():

    def __new__(cls, fromrational=None, **kwargs):
        obj = super().__new__(cls)
        if 'prime' in kwargs:
            obj.prime = kwargs['prime']
        else:
            obj.prime = 2
        if 'significand' in kwargs:
            if fromrational is not None:
                raise ValueError("Both rational value and significand given, "
                                 "expected only one or none.")
            if 'exponent' not in kwargs:
                raise ValueError('Significand given, but exponent not.')
            obj.significand = kwargs['significand']
            obj.exponent = kwargs['exponent']
            return obj
        elif fromrational is not None:
            if isinstance(fromrational, int):
                obj.exponent = obj._valuationFromInt(fromrational)
                obj.significand = (fromrational // (obj.prime**obj.exponent)) \
                    % (obj.prime ** _PADIC_PRECISION)
                return obj
            if isinstance(fromrational, Fraction):
                num = fromrational.numerator
                dem = fromrational.denominator
                numval = obj._valuationFromInt(num)
                if numval == 0:  # negative valuation
                    demval = obj._valuationFromInt(dem)
                    unitdem = dem // (obj.prime ** demval)
                    obj.exponent = -demval
                    unitfactor = num
                else:
                    demval = 0  # we use that fractions are reduced
                    obj.exponent = numval
                    unitfactor = num // (obj.prime ** numval)
                    unitdem = dem
                deminverse = _invmod(unitdem, obj.prime, _PADIC_PRECISION)
                obj.significand = (unitfactor*deminverse) \
                    % (obj.prime ** _PADIC_PRECISION)
                return obj
            else:
                t = type(fromrational)
                raise TypeError(f"Cannot create p-adic from type {t}")
        else:
            obj.significand = 0
            obj.exponent = _MAX_PADIC_EXPONENT
            return obj

    @classmethod
    def NaN(cls, prime=2):
        return PAdicFloat(exponent=_MIN_PADIC_EXPONENT-1,
                          significand=0,
                          prime=prime)

    @classmethod
    def inf(cls, prime=2):
        return PAdicFloat(exponent=_MIN_PADIC_EXPONENT-1,
                          significand=1,
                          prime=prime)

    def _valuationFromInt(self, i):
        val = 0
        if i == 0:
            return _MAX_PADIC_EXPONENT
        while i % self.prime == 0:
            i = i//self.prime
            val += 1
        return val

    def _significantDisplay(self):
        coeffs = [0] * _DISPLAY_CAP
        n = self.significand
        for j in _DISPLAY_RANGE:
            n, r = divmod(n, self.prime)
            coeffs[j] = r
            if n == 0:
                break
        return tuple(coeffs)

    def __repr__(self):
        return f"({self.exponent}, {self.significand})_{self.prime}"

    def __str__(self):
        if self.iszero():
            return "0"
        if self.isNaN():
            return "NaN"
        if self.isinf():
            return "∞" if _USE_UNICODE else "inf"
        digitStr = ' '.join([str(d)
                             for d in self._significantDisplay()])
        return f"({digitStr}){self.prime}**{self.exponent}"

    def iszero(self):
        return self.exponent == _MAX_PADIC_EXPONENT and self.significand == 0

    def isinf(self):
        return self.exponent == _MIN_PADIC_EXPONENT - 1 and \
            self.significand == 1

    def isNaN(self):
        return self.exponent == _MIN_PADIC_EXPONENT - 1 and \
            self.significand == 0

    def _report(self):
        print(f"z: {self}")
        print(f"z is zero? {self.iszero()}")
        print(f"z is inf? {self.isinf()}")
        print(f"z is NaN? {self.isNaN()}")

    def normalize(self):
        if self.exponent < _MIN_PADIC_EXPONENT and self.significand == 0:
            exponent = _MIN_PADIC_EXPONENT - 1
            significand = 0  # NaN
        elif self.exponent < _MIN_PADIC_EXPONENT:
            exponent = _MIN_PADIC_EXPONENT - 1
            significand = 1  # inf
        elif self.exponent > _MAX_PADIC_EXPONENT:
            exponent = _MAX_PADIC_EXPONENT
            significand = 0  # zero
        elif self.significand == 0:
            exponent = _MAX_PADIC_EXPONENT
            significand = 0  # zero
        else:
            exponent = self._valuationFromInt(self.significand)
            significand = (self.significand // (self.prime ** exponent)) \
                % (self.prime ** _PADIC_PRECISION)
            exponent += self.exponent
        return PAdicFloat(exponent=exponent,
                          significand=significand,
                          prime=self.prime)

    def __eq__(self, other):
        if not isinstance(other, PAdicFloat):
            try:
                other = PAdicFloat(other, prime=self.prime)
            except TypeError:
                return NotImplemented
        a = self.normalize()
        b = other.normalize()
        if a.isNaN() or b.isNaN():
            return False
        return (a.exponent == b.exponent) and \
               (a.significand == b.significand) and \
               (a.prime == b.prime)

    def __bool__(self):
        return not self.iszero()

    def __mul__(self, other):
        if not isinstance(other, PAdicFloat):
            try:
                other = PAdicFloat(other, prime=self.prime)  # tries to convert
            except TypeError:
                return NotImplemented
        if (self.prime != other.prime):
            return NotImplemented
        a = self.normalize()
        b = other.normalize()
        if a.isNaN() or b.isNaN():
            return PAdicFloat.NaN(prime=self.prime)
        if a.isinf() or b.isinf():
            if a.iszero() or b.iszero():
                return PAdicFloat.NaN(prime=self.prime)
            return PAdicFloat.inf(prime=self.prime)
        newexp = a.exponent + b.exponent
        if newexp > _MAX_PADIC_EXPONENT:  # underflow
            newsigs = 0
            newexp = _MAX_PADIC_EXPONENT
        else:
            newsigs = (a.significand*b.significand) \
                % (self.prime ** _PADIC_PRECISION)
        return PAdicFloat(exponent=newexp,
                          significand=newsigs,
                          prime=self.prime)

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        if not isinstance(other, PAdicFloat):
            try:
                other = PAdicFloat(other, prime=self.prime)  # tries to convert
            except TypeError:
                return NotImplemented
        if (self.prime != other.prime):
            return NotImplemented
        a = self.normalize()
        b = other.normalize()
        if a.isNaN() or b.isNaN():
            return PAdicFloat.NaN(prime=self.prime)
        if b.iszero():
            if a.iszero():
                return PAdicFloat.NaN(prime=self.prime)
            return PAdicFloat.inf(prime=self.prime)
        if b.isinf():
            if a.isinf():
                return PAdicFloat.NaN(prime=self.prime)
            return PAdicFloat(0, prime=self.prime)
        newexp = a.exponent - b.exponent
        if newexp < _MIN_PADIC_EXPONENT:  # overflow
            newsigs = 1
            newexp = _MIN_PADIC_EXPONENT - 1
        else:
            c = _invmod(b.significand, self.prime, _PADIC_PRECISION)
            newsigs = (a.significand * c) % (self.prime ** _PADIC_PRECISION)
        return PAdicFloat(exponent=newexp,
                          significand=newsigs,
                          prime=self.prime)

    def __rtruediv__(self, other):
        return PAdicFloat(other, prime=self.prime) / self

    def __add__(self, other):
        if not isinstance(other, PAdicFloat):
            try:
                other = PAdicFloat(other, prime=self.prime)  # tries to convert
            except TypeError:
                return NotImplemented
        if (self.prime != other.prime):
            return NotImplemented
        a = self.normalize()
        b = other.normalize()
        if a.isNaN() or b.isNaN():
            return PAdicFloat.NaN(prime=self.prime)
        if a.isinf() or b.isinf():
            return PAdicFloat.inf(prime=self.prime)
        if (a.exponent == b.exponent):
            # cancellation may happen
            v = self._valuationFromInt(a.significand + b.significand)
            if v > _MAX_PADIC_EXPONENT - a.exponent:  # overflow
                newexp = _MAX_PADIC_EXPONENT
                newsigs = 0
            else:
                newexp = a.exponent + v
                newsigs = (a.significand + b.significand) // (self.prime ** v)
        elif a.exponent > b.exponent:
            # non-archimedianess kicks in!
            newexp = b.exponent
            newsigs = (a.significand
                       // (self.prime**(a.exponent - b.exponent))) \
                + b.significand
        else:
            # non-archimedianess kicks in!
            newexp = a.exponent
            newsigs = (b.significand
                       // (self.prime**(b.exponent - a.exponent))) \
                + a.significand
        return PAdicFloat(significand=newsigs,
                          exponent=newexp,
                          prime=self.prime)

    def __radd__(self, other):
        return self+other

    def __neg__(self):
        return self*-1

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -(self - other)

    def __invert__(self):
        return 1/self
