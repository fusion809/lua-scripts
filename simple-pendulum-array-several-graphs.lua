--
--------------------------------------------------------------------------------
--         File:  simple-pendulum.lua
--
--        Usage:  ./simple-pendulum.lua
--
--  Description:  Solves the problem of the double pendulum using Runge-Kutta 4th order method
--                minth is 7.344e-12 if N = 1e8 and tN = 1.5 (took 114.26s).
--      Options:  ---
-- Requirements:  ---
--         Bugs:  ---
--        Notes:  ---
--       Author:  Brenton Horne (), <brentonhorne77@gmail.com>
-- Organization:  
--      Version:  1.0
--      Created:  16/10/17
--     Revision:  ---
--------------------------------------------------------------------------------
-- Define parameters
local t0         = 0                            -- Initial time
local tN         = 4                            -- Finishing time
local theta0     = 0                            -- Angle from positive x axis (theta) at t0
local dtheta0    = 0                            -- Change rate in theta at t0
local N          = 1e7                          -- Number of steps
local g          = 9.8                          -- Acceleration due to gravity
local l          = 1                            -- Length of pendulum

-- Initiate problem variables
local t          = {}                           -- Initialize t array
local theta      = {}                           -- Initialize theta array
local dtheta     = {}                           -- Initialize dtheta array
local d2theta    = {}                           -- Initialize d2theta array
local d3theta    = {}                           -- Initialize d3theta array
local d4theta    = {}                           -- Initialize d4theta array
local d5theta    = {}                           -- Initialize d5theta array
local d6theta    = {}                           -- Initialize d6theta array
t[1]             = t0                           -- Initiate t[1] variable
theta[1]         = theta0                       -- Initiate theta[1] variable
dtheta[1]        = dtheta0                      -- Initiate dtheta[1] variable
d2theta[1]       = - g / l * math.cos(theta[1]) -- Initiate d2theta[1] variable
d3theta[1]       = g / l * dtheta[1] * math.sin(theta[1]) -- Initiate d3theta[1] variable
d4theta[1]       = - 3 * ( g^2 ) / ( l^2 ) * math.cos(theta[1]) * math.sin(theta[1]) -- Initiate 4theta[1] variable
d5theta[1]       = 3 * ( g^2 ) / ( l^2 ) * dtheta[1] * ((math.sin(theta[1]))^2 - (math.cos(theta[1]))^2)
d6theta[1]       = 3 * ( g^3 ) / ( l^3 ) * ( ( math.cos(theta[1]) ) * ( math.cos(2*theta[1]) ) - 4 * ( math.sin(theta[1]) ) * ( math.sin(2 * theta[1]) ) )
h                = (tN - t0) / N                -- define step size
minth            = theta0 + math.pi             -- Define minth at t0

-- Define the d2theta/dt2 = f(g, l, t, theta, dtheta) function
-- Essentially the RHS of the problem
function f(g, l, t, theta, dtheta)
    K = - (g/l) * math.cos(theta)
    return K
end

local k1, l1, k2, l2, k3, l3, k4, l4
-- Loop over time
for i = 1,N do
    -- First approximation 
    k1           = h * f(g, l, t[i],       theta[i],            dtheta[i])
    l1           = h * dtheta[i];

    -- Second approximation
    k2           = h * f(g, l, t[i] + h/2, theta[i] + 1/2 * l1, dtheta[i] + 1/2 * k1)
    l2           = h * (dtheta[i] + 1/2 * k1)

    -- Third approximation
    k3           = h * f(g, l, t[i] + h/2, theta[i] + 1/2 * l2, dtheta[i] + 1/2 * k2)
    l3           = h * (dtheta[i] + 1/2 * k2)

    -- Fourth approximation
    k4           = h * f(g, l, t[i] + h,   theta[i] + l3,       dtheta[i] + k3)
    l4           = h * (dtheta[i] + 1/2 * k3)

    -- Updating variables
    t[i+1]       = t[i] + h
    dtheta[i+1]  = dtheta[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    theta[i+1]   = theta[i]  + (1 / 6) * (l1 + 2 * l2 + 2 * l3 + l4)

    d2theta[i+1] = - g / l * math.cos(theta[i+1])
    d3theta[i+1] = g / l * dtheta[i+1] * math.sin(theta[i+1])
    d4theta[i+1] = - 3 * ( g^2 ) / ( l^2 ) * math.cos(theta[i+1]) * math.sin(theta[i+1])
    d5theta[i+1] = 3 * ( g^2 ) / ( l^2 ) * dtheta[i+1] * ((math.sin(theta[i+1]))^2 - (math.cos(theta[i+1]))^2)
    d6theta[i+1] = 3 * ( g^3 ) / ( l^3 ) * ( ( math.cos(theta[i+1]) ) * ( math.cos(2*theta[i+1]) ) - 4 * ( math.sin(theta[i+1]) ) * ( math.sin(2 * theta[i+1]) ) )
    -- Determining pi + theta; at minima it is 0
    diff         = math.abs(theta[i] + math.pi)
 
    -- Checking if diff is smaller than the smallest minth so far and updating it if it is 
    if (diff < minth) then
         minth   = diff
    end
         
end

print(minth)

--Plot theta against t
gp = require("gnuplot")           -- Uses this library https://bitbucket.org/lucashnegri/lua-gnuplot/raw/34136a285a8820f31e1c63b7c81aa1b70a4b60ec/gnuplot.lua

--[[local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data = {
        gp.array {
            {
                t,
                theta
            },
        },
    }
}:plot('01-simple-pendulum-theta-by-t.png')

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                dtheta
            },
        },
    }
}:plot('02-simple-pendulum-dtheta-by-t.png')

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                d2theta
            },
        },
    }
}:plot('03-simple-pendulum-d2theta-by-t.png')

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                d3theta
            },
        },
    }
}:plot('04-simple-pendulum-d3theta-by-t.png')

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                d4theta
            },
        },
    }
}:plot('05-simple-pendulum-d4theta-by-t.png')

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                d5theta
            },
        },
    }
}:plot('06-simple-pendulum-d5theta-by-t.png')]]--

local g = gp{
    width        = 1600,
    height       = 900,
    xlabel       = "X axis",
    ylabel       = "Y axis",
    key          = "top left",
    consts       = {
        gamma    = 2.5
    },

    data         = {
        gp.array {
            {
                t,
                d6theta
            },
        },
    }
}:plot('07-simple-pendulum-d6theta-by-t.png')


--[[Plot theta against t
lp = require("lua-plot")
p = lp.plot{}
p:AddSeries(t, d2theta)
p:Show()
io.read()]]--
