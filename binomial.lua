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
--
local n
io.write("Enter an integer value for n ")
io.flush()
n = io.read()

local k
io.write("Enter an integer value for k ")
io.flush()
k = io.read()

local C = n

for i = 1,k-1 do
    C = (n-i)*C / (i+1)
end

print(C)
