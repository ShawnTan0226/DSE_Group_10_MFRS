import math
# Input parameters
taper_ratio = float(input("Enter the taper ratio: "))
aspect_ratio = float(input("Enter the aspect ratio: "))
wing_sweep = math.radians(float(input("Enter the wing sweep (in degrees): ")))
wingspan = float(input("Enter the wingspan: "))



# Calculate mean aerodynamic chord
mean_aero_chord = (2/3) * wingspan * ((1 + taper_ratio + taper_ratio**2) / (1 + taper_ratio)) / aspect_ratio

# Correct for wing sweep
mean_aero_chord /= math.cos(wing_sweep)

# Output result
print(f"The mean aerodynamic chord is {mean_aero_chord:.2f} units.")

#Find x_mac

