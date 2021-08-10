from fractions import Fraction

_PADIC_PRECISION = 65  # we're storing 65 p-adic digits
_MAX_PADIC_EXPONENT = 64  # max valuation is 64, min is -63
_MIN_PADIC_EXPONENT = 1 - _MAX_PADIC_EXPONENT
_DISPLAY_CAP = 10  # when printed, show first 10 digits
_USE_UNICODE = True  # should we print ∞ as unicode on screen?

_PRECISION_RANGE = range(0, _PADIC_PRECISION)
_DISPLAY_RANGE = range(0, _DISPLAY_CAP)


def _invmod(a, p, power):  # is there a faster way to do this?
    # Hensel's lemma
    c = pow(a, -1, p)
    b = c
    for i in range(2, power+1):
        b = (b - (a*b - 1)*c) % p**i
    return b


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
                obj.exponent = cls._valuationFromInt(fromrational, obj.prime)
                obj.significand = (fromrational // (obj.prime**obj.exponent)) \
                    % (obj.prime ** _PADIC_PRECISION)
                return obj
            if isinstance(fromrational, Fraction):
                num = fromrational.numerator
                dem = fromrational.denominator
                numval = cls._valuationFromInt(num, obj.prime)
                if numval == 0:  # negative valuation
                    demval = cls._valuationFromInt(dem, obj.prime)
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

    @classmethod
    def _valuationFromInt(cls, i, p):
        val = 0
        if i == 0:
            return _MAX_PADIC_EXPONENT
        while i % p == 0:
            i = i // p
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
            exponent = self.__class__._valuationFromInt(self.significand,
                                                        self.prime)
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
        return (a-b).iszero()

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
            v = self.__class__._valuationFromInt(a.significand + b.significand,
                                                 self.prime)
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
                       * (self.prime**(a.exponent - b.exponent))) \
                + b.significand
        else:
            # non-archimedianess kicks in!
            newexp = a.exponent
            newsigs = (b.significand
                       * (self.prime**(b.exponent - a.exponent))) \
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


def plog(padicNumber):
    t = (padicNumber - 1).normalize()
    texp = t.exponent
    tclean = t / (t.prime ** texp)
    if texp < 1:
        raise ValueError("p-adic logarithm only defined for 1 + pZp.")
    log = PAdicFloat(0, prime=t.prime)
    summand = t
    powersig = tclean
    powerexp = texp
    target = (_MAX_PADIC_EXPONENT // texp) + 1
    for i in range(2, target+1):
        log = (log + summand).normalize()
        powersig *= tclean
        powerexp += texp
        iexp = PAdicFloat._valuationFromInt(i, t.prime)
        ifactor = i // (t.prime ** iexp)
        factor = PAdicFloat(Fraction(pow(-1, i-1), ifactor), prime=t.prime)
        scaling = PAdicFloat(t.prime ** (powerexp - iexp), prime=t.prime)
        summand = ((powersig * factor) * scaling).normalize()
    return log


def pexp(padicNumber):
    t = padicNumber.normalize()
    texp = t.exponent
    tclean = t / (t.prime**texp)
    if t.prime == 2 and texp < 2:
        raise ValueError("p-adic exponential only defined for 4Z_2")
    if t.prime != 2 and texp < 1:
        raise ValueError("p-adic exponential only defined for pZp")
    exp = PAdicFloat(0, prime=t.prime)
    summand = PAdicFloat(1, prime=t.prime)
    powersig = PAdicFloat(1, prime=t.prime)
    powerexp = 0
    # to determine bounts, we use that val_p(n!) ~ n/(p-1)
    # thus, n*val_p(t) - val_p(n!) ~ n*(val_p(t) - 1/(p-1))
    # assymptotically.
    # hence, we sum only up to n = (precision / (vap_p(t) - 1/(p-1))) + 1
    target = int(_MAX_PADIC_EXPONENT / (t.exponent - 1/(t.prime - 1))) + 1
    factexp = 0
    factsig = PAdicFloat(1, prime=t.prime)
    for i in range(1, target+1):
        exp = (exp + summand).normalize()
        powerexp += texp
        powersig *= tclean
        iexp = PAdicFloat._valuationFromInt(i, t.prime)
        factexp += iexp
        factsig /= (i // (t.prime**iexp))
        scaling = PAdicFloat(t.prime**(powerexp-factexp),
                             prime=t.prime).normalize()
        summand = (scaling * (powersig * factsig)).normalize()
    return exp
