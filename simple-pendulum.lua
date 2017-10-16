--
--------------------------------------------------------------------------------
--         File:  simple-pendulum.lua
--
--        Usage:  ./simple-pendulum.lua
--
--  Description:  Solves the problem of the double pendulum using Runge-Kutta 4th order method
--
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
local t0      = 0
local tN      = 100
local theta0  = 0
local dtheta0 = 0
local N       = 1e9
local g       = 9.8
local l       = 1

-- Initiate problem variables
local t       = t0
local theta   = theta0
local dtheta  = dtheta0
local h       = (tN - t0) / N
local minth   = math.pi

-- Define the d2theta/dt2 = f(g, l, t, theta, dtheta) function
-- Essentially the RHS of the problem
function f(g, l, t, theta, dtheta)
    K = - (g/l) * math.cos(theta)
    return K
end

-- Loop over time
for i = 1,N do
    -- First approximation
    k1     = h * f(g, l, t, theta, dtheta)
    l1     = h * dtheta;

    -- Second approximation
    k2     = h * f(g, l, t + h/2, theta + 1/2 * l1, dtheta + 1/2 * k1)
    l2     = h * (dtheta + 1/2 * k1)

    -- Third approximation
    k3     = h * f(g, l, t + h/2, theta + 1/2 * l2, dtheta + 1/2 * k2)
    l3     = h * (dtheta + 1/2 * k2)

    -- Fourth approximation
    k4     = h * f(g, l, t + h, theta + l3, dtheta + k3)
    l4     = h * (dtheta + 1/2 * k3)

    -- Updating variables
    t      = t + h
    dtheta = dtheta + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    theta  = theta + (1 / 6) * (l1 + 2 * l2 + 2 * l3 + l4)

    -- Determining pi + theta; at minima it is 0
    diff   = math.abs(theta + math.pi)

    -- Checking if diff is smaller than the smallest minth so far and updating it if it is 
    if (diff < minth) then
         minth = diff
    end
         
end

-- Print it
print(minth)
