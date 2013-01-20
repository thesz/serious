-- |MAtrix.hs
--
-- Defines sparse matrix in hypersparse format (see the work of Buluc) and provides some operations
-- with them.
--
-- Copyright (C) 2013 Serguey Zefirov

module Math.Sparse.Matrix(
          SM
        , empty
        , singleton
        , (.+), (.-), (.*)
        , spMMV
        ) where

import Data.List (sortBy)

-- |A class to check for zero values and to provide them.
-- Useful to define tests for floating point numbers. We also provide obvious one for integers.
class Zero a where
        -- |Is value a zero?
        isZero :: a -> Bool
        -- |Zero value. isZero zeroValue == True
        zeroValue :: a

instance Zero Int where
        isZero = (==0)
        zeroValue = 0

-- |Indexed element.
type IElem a = (Int, Int, a)

-- |And a list of indexed elements.
type IElems a = [IElem a]

-- |The hypersparse matrix, an abstract type.
data SM a = SM {
        -- |Sorted list of elements.
          smElements :: IElems a
        }
        deriving Show

-- |Empty matrix.
empty :: SM a
empty = SM []

-- |Singleton matrix.
singleton :: Int -> Int -> a -> SM a
singleton row col a = SM [(row, col, a)]

-- |Helper - merges elements of two matrices with the function provided.
-- Zero elements are filtered out.
merge :: (Num a, Zero a) => (a -> a -> a) -> IElems a -> IElems a -> IElems a
merge f [] bs = bs
merge f as [] = as
merge f as@(x@(ia,ja,a):ijas) bs@(y@(ib,jb,b):ijbs) = case compare (ia,ja) (ib, jb) of
        LT -> x : merge f ijas bs
        GT -> y : merge f as ijbs
        EQ
                | isZero r -> merge f ijas ijbs
                | otherwise -> (ia,ja,r) : merge f ijas ijbs
        where
                r = f a b

-- |The addition.
(.+) :: (Num a, Zero a) => SM a -> SM a -> SM a
infixl 6 .+
SM ea .+ SM eb = SM $ merge (+) ea eb

-- |The subtraction.
(.-) :: (Num a, Zero a) => SM a -> SM a -> SM a
infixl 6 .-
SM ea .- SM eb = SM $ merge (-) ea eb

sortByIndices :: IElems a -> IElems a
sortByIndices = sortBy compareIndices
        where
                compareIndices (i1,j1,_) (i2,j2,_) = compare (i1,j1) (i2,j2)

-- |Helper - transposes elements.
transposeElems :: IElems a -> IElems a
transposeElems es = sortByIndices $ map (\(i,j,a) -> (j,i,a)) es

-- |Transposition.
transpose :: SM a -> SM a
transpose (SM es) = SM $ transposeElems es

-- |The multiplication.
-- This is relatively tricky.
-- (ia,k,a) should be multiplied to all (k,jb,b) and summed, if ia == jb.
(.*) :: (Num a, Zero a) => SM a -> SM a -> SM a
infixl 7 .*
SM ea .* SM eb = SM $ reduce $ sortByIndices outer
        where
                reduce (a1@(i1,j1,x1) : rest@((i2,j2,x2) : ijxs))
                        | i1 == i2 && j1 == j2 =
                                if isZero y then reduce ijxs else reduce ((i1,j1,y):ijxs)
                        | otherwise = a1 : reduce rest
                        where
                                y = x1+x2
                reduce ijxs = ijxs
                outer = [(ia,jb,a*b) | (ia,ja,a) <- ea, (ib,jb,b) <- eb, ja == ib]

-- |Sparse vector multiplication.
spMMV :: (Num a, Zero a) => SM a -> [(Int, a)] -> [(Int, a)]
spMMV sm sv = shrink r
        where
                shrink = map (\(i,_,a) -> (i,a))
                expand = SM . sortByIndices . map (\(i,a) -> (i,1,a))
                (SM r) = sm .* expand sv
