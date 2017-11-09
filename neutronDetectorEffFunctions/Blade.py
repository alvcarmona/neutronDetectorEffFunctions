

class Blade:
    backscatter = 0.0
    transmission = 0.0
    substrate = 0.0
    inclination = 0

    def __init__(self, back,tra, sub, incl, parent=None):
        self.backscatter = back
        self.transmission = tra
        self.substrate = sub
        self. inclination = incl


