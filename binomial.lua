--
--------------------------------------------------------------------------------
--         File:  binomial.lua
--
--        Usage:  ./binomial.lua
--
--  Description:  Calculates the binomial coefficient
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
-- Determines n choose k

-- Get n
local n
io.write("Enter an integer value for n ")
io.flush()
n = io.read()

-- Get k
local k
io.write("Enter an integer value for k ")
io.flush()
k = io.read()

-- Initiate C
local C = n

-- Loop over i
for i = 1,k-1 do
    -- Calculate C for ith iteration
    C = (n-i)*C / (i+1)
end

-- C now equals the binomial coefficient
print(C)
