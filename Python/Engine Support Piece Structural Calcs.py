in_to_m = 0.0254

m_engine = 10.611 # kg
g = 9.81 # m/s^2

weight_engine = m_engine * g # N
length = 8.0 * in_to_m # m
depth = 0.75 * in_to_m # m

normalStress = weight_engine/(length * depth) / 10**6 # MPa
yieldStress = 30 # MPa

if normalStress > yieldStress:
    print('Part will yield') 
    
print(f"Normal stress: {normalStress} MPa, Yield Stress: {yieldStress} MPa")