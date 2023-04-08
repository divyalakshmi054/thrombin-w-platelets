# packages -
using BSTModelKit
using CSV
using DataFrames
using DifferentialEquations
using GlobalSensitivity
using NumericalIntegration

# paths -
_PATH_TO_MODEL = joinpath(pwd(),"model")
_PATH_TO_DATA = joinpath(pwd(),"data")
_PATH_TO_TMP = joinpath(pwd(),"tmp")