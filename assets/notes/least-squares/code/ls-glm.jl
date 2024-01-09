# This file was generated, do not modify it. # hide
using GLM
ols = lm(@formula(MPG ~ Disp + HP + WT), mtcars)