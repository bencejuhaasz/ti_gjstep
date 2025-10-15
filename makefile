# ----------------------------
# Makefile Options
# ----------------------------

NAME = GJSTEP
ICON = icon.png
DESCRIPTION = "GAUSS-JORDAN SOLVER"
COMPRESSED = NO

CFLAGS = -Wall -Wextra -Oz
CXXFLAGS = -Wall -Wextra -Oz

# ----------------------------

include $(shell cedev-config --makefile)
