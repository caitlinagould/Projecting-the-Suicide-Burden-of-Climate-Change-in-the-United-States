#!/bin/bash

CONFIG="configuration.yaml"
DEGREE="D6" # D1 or D2 or D3 or D4 or D5 or D6
POINT_MODE="FALSE"
POPULATION_YEAR="FUTURE" # PRESENT or FUTURE

INCOME_YEAR="PRESENT" # PRESENT or FUTURE
DISCOUNT_RATE=3
DISCOUNT_YEAR=2015


Rscript simulate/estimate.R  ${CONFIG} ${DEGREE} ${POINT_MODE} ${POPULATION_YEAR} ${INCOME_YEAR} ${DISCOUNT_RATE} ${DISCOUNT_YEAR}
Rscript simulate/testConvergence.R  ${CONFIG} ${DEGREE} ${POINT_MODE} ${POPULATION_YEAR} ${INCOME_YEAR} ${DISCOUNT_RATE} ${DISCOUNT_YEAR}

