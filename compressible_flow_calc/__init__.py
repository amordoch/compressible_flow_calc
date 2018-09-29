from compressible_flow_calc import iteration
from compressible_flow_calc import calc

# Import pint if it's installed. If not, create a dummy class which returns 1
# for any requested unit.
try:
    import pint
except ModuleNotFoundError:
    class UnitDisabler:
        def __getattr__(self, item):
            return 1

    unit = UnitDisabler()

    def Q(x, y):
        return x
else:
    # Instantiate a global unit registry
    unit = pint.UnitRegistry()
    Q = unit.Quantity
    pint.set_application_registry(unit)

__CFCVERSION__ = '0.1.3'