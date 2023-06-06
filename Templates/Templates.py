#####################################     TEMPLATES     ###############################


# Template for graphs:
fig, ax = plt.subplots()

# Possible colours: 1) blue, 2) red, 3) green, 4) yellow, 5)
ax.plot(x, y, linestyle='-', label="Label", color='SPECIFY COLOUR!')

ax.set_xlabel('Axis Label')
ax.set_ylabel('Axis Label')
ax.minorticks_on()
ax.grid(which='major', color='grey', linestyle='-')
ax.grid(which='minor', color='lightgrey', linestyle='dotted')
ax.axhline(y=0, xmin=0, xmax=1, linewidth=2, color='black')
# plt.gca().invert_yaxis()
ax.legend()
plt.title("Title")
plt.show()