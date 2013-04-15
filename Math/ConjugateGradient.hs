---------------------------------------------------------------------------------
-- |
-- Module      :  Math.ConjugateGradient
-- Copyright   :  (c) Levent Erkok
-- License     :  BSD3
-- Maintainer  :  erkokl@gmail.com
-- Stability   :  stable
--
-- (The linear equation solver library is hosted at <http://github.com/LeventErkok/conjugateGradient>.
-- Comments, bug reports, and patches are always welcome.)
--
-- Sparse matrix linear-equation solver, using the conjugate gradient algorithm. Note that the technique only
-- applies to matrices that are:
--
--              * Symmetric
--
--              * Positive-definite
--
--  See <http://en.wikipedia.org/wiki/Conjugate_gradient_method> for details.
--
--  The conjugate gradient method can handle very large sparse matrices, where direct
--  methods (such as LU decomposition) are way too expensive to be useful in practice.
--  Such large sparse matrices arise naturally in many engineering problems, such as
--  in ASIC placement algorithms and when solving partial differential equations.
--
-- Here's an example usage, for the simple system:
--
-- @
--       4x +  y = 1
--        x + 3y = 2
-- @
--
-- >>> import Data.IntMap
-- >>> let a = fromList [(0, fromList [(0, 4), (1, 1)]), (1, fromList [(0, 1), (1, 3)])] :: SM Double
-- >>> let b = fromList [(0, 1), (1, 2)] :: SV Double
-- >>> let g = mkStdGen 12345
-- >>> let (_, x) = solveCG g 2 a b
-- >>> putStrLn $ showSolution 4 2 a b x
--       A       |   x    =   b   
-- --------------+----------------
-- 4.0000 1.0000 | 0.0909 = 1.0000
-- 1.0000 3.0000 | 0.6364 = 2.0000
---------------------------------------------------------------------------------

module Math.ConjugateGradient(
          -- * Types
          -- $typeInfo
            SV, SM
          -- * Sparse operations
          , sMulSV, sMulSM, addSV, subSV, dotSV, mulSMV, normSV
          -- * Conjugate-Gradient solver
          , solveCG
          -- * Displaying solutions
          , showSolution
        ) where

import Data.List                   (intercalate)
import Data.Maybe                  (fromMaybe)
import qualified Data.IntMap as IM (IntMap, lookup, map, unionWith, intersectionWith, fold, fromList)
import System.Random
import Numeric

-- | A sparse vector containing elements of type 'a'. (For our purposes, the elements will be either 'Float's or 'Double's.)
type SV a = IM.IntMap a

-- | A sparse matrix is an int-map containing sparse row-vectors. Again, only put in rows that have a non-@0@ element in them for efficiency.
type SM a = IM.IntMap (SV a)

---------------------------------------------------------------------------------
-- Sparse vector/matrix operations
---------------------------------------------------------------------------------

-- | Look-up a value in a sparse-vector.
vLookup :: Num a => SV a -> Int -> a
vLookup m k = fromMaybe 0 (k `IM.lookup` m)

-- | Look-up a value in a sparse-matrix.
mLookup :: Num a => SM a -> (Int, Int) -> a
mLookup m (i, j) = maybe 0 (`vLookup` j) (i `IM.lookup` m)

-- | Multiply a sparse-vector by a scalar.
sMulSV :: RealFloat a => a -> SV a -> SV a
sMulSV s = IM.map (s *)

-- | Multiply a sparse-matrix by a scalar.
sMulSM :: RealFloat a => a -> SM a -> SM a
sMulSM s = IM.map (s `sMulSV`)

-- | Add two sparse vectors.
addSV :: RealFloat a => SV a -> SV a -> SV a
addSV = IM.unionWith (+)

-- | Subtract two sparse vectors.
subSV :: RealFloat a => SV a -> SV a -> SV a
subSV v1 v2 = addSV v1 (IM.map ((-1)*) v2)

-- | Dot product of two sparse vectors.
dotSV :: RealFloat a => SV a -> SV a -> a
dotSV v1 v2 = IM.fold (+) 0 $ IM.intersectionWith (*) v1 v2

-- | Multiply a sparse matrix (nxn) with a sparse vector (nx1), obtaining a sparse vector (nx1).
mulSMV :: RealFloat a => SM a -> SV a -> SV a
mulSMV m v = IM.map (`dotSV` v) m

-- | Norm of a sparse vector. (Square-root of its dot-product with itself.)
normSV :: RealFloat a => SV a -> a
normSV = sqrt . IM.fold (\e s -> e*e + s) 0

-- | Conjugate Gradient Solver for the system @Ax=b@. See: <http://en.wikipedia.org/wiki/Conjugate_gradient_method>.
--
-- NB. Assumptions on the input:
--
--    * The @A@ matrix is symmetric and positive definite.
--
--    * The indices start from @0@ and go consecutively up-to @n-1@. (Only non-@0@ value/row
--      indices has to be present, of course.)
--
-- For efficiency reasons, we do not check for either property. (If these assumptions are
-- violated, the algorithm will still produce a result, but not the one you expected!)
--
-- We perform either @10^6@ iterations of the Conjugate-Gradient algorithm, or until the error
-- factor is less than @1e-10@. The error factor is defined as the difference of the norm of
-- the current solution from the last one, as we go through the iteration. See
-- <http://en.wikipedia.org/wiki/Conjugate_gradient_method#Convergence_properties_of_the_conjugate_gradient_method>
-- for a discussion on the convergence properties of this algorithm.
solveCG :: (RandomGen g, RealFloat a, Random a)
        => g          -- ^ The seed for the random-number generator.
        -> Int        -- ^ Number of variables.
        -> SM a       -- ^ The @A@ sparse matrix (@nxn@).
        -> SV a       -- ^ The @b@ sparse vector (@nx1@).
        -> (a, SV a)  -- ^ The final error factor, and the @x@ sparse matrix (@nx1@), such that @Ax = b@.
solveCG g n a b = cg a b x0
  where rs = take n (randomRs (0, 1) g)
        x0 = IM.fromList [p | p@(_, j) <- zip [0..] rs, j /= 0]

-- | The Conjugate-gradient algorithm. Our implementation closely follows the
-- one given here: <http://en.wikipedia.org/wiki/Conjugate_gradient_method#Example_code_in_Matlab>
cg :: RealFloat a => SM a -> SV a -> SV a -> (a, SV a)
cg a b x0 = cgIter (1000000 :: Int) (norm r0) r0 r0 x0
 where r0 = b `subSV` (a `mulSMV` x0)
       cgIter 0 eps _ _ x = (eps, x)
       cgIter i eps p r x
        -- Stop if the square of the error is less than 1e-20, i.e.,
        -- if the error itself is less than 1e-10.
        | eps' < 1e-20 = (eps', x')
        | True         = cgIter (i-1) eps' p' r' x'
        where ap    = a `mulSMV` p
              alpha = eps / ap `dotSV` p
              x'    = x `addSV` (alpha `sMulSV` p)
              r'    = r `subSV` (alpha `sMulSV` ap)
              eps'  = norm r'
              p'    = r' `addSV` ((eps' / eps) `sMulSV` p)
       norm = IM.fold (\e s -> e*e + s) 0  -- square of normSV, but no need for expensive square-root

-- | Display a solution in a human-readable form. Needless to say, only use this
-- method when the system is small enough to fit nicely on the screen.
showSolution :: RealFloat a
             => Int   -- ^ Precision: Use this many digits after the decimal point.
             -> Int   -- ^ Number of variables, @n@
             -> SM a  -- ^ The @A@ matrix, @nxn@
             -> SV a  -- ^ The @b@ matrix, @nx1@
             -> SV a  -- ^ The @x@ matrix, @nx1@, as returned by 'solveCG', for instance.
             -> String
showSolution prec n ma vb vx = intercalate "\n" $ header ++ res
  where res   = zipWith3 row a x b
        range = [0..n-1]
        sf d = showFFloat (Just prec) d ""
        a = [[sf (ma `mLookup` (i, j)) | j <- range] | i <- range]
        x = [sf (vx `vLookup` i) | i <- range]
        b = [sf (vb `vLookup` i) | i <- range]
        cellWidth = maximum (0 : map length (concat a ++ x ++ b))
        row as xv bv = unwords (map pad as) ++ " | " ++ pad xv ++ " = " ++ pad bv
        pad s  = reverse $ take (length s `max` cellWidth) $ reverse s ++ repeat ' '
        center l s = let extra         = l - length s
                         (left, right) = (extra `div` 2, extra - left)
                     in  replicate left ' ' ++ s ++ replicate right ' '
        header = case res of
                   []    -> ["Empty matrix"]
                   (r:_) -> let l = length (takeWhile (/= '|') r)
                                h =  center (l-1) "A"  ++ " | "
                                  ++ center cellWidth "x" ++ " = " ++ center cellWidth "b"
                                s = replicate l '-' ++ "+" ++ replicate (length r - l - 1) '-'
                            in [h, s]

{- $typeInfo
We represent sparse matrices and vectors using 'IM.IntMap's. In a sparse vector, we only populate those elements that are non-@0@.
In a sparse matrix, we only populate those rows that contain a non-@0@ element. This leads to an efficient representation for
sparse matrices and vectors, where the space usage is proportional to number of non-@0@ elements. Strictly speaking, putting non-@0@ elements
would not break the algorithms we use, but clearly they would be less efficient.

Indexings starts at @0@, and is assumed to be non-negative, corresponding to the row numbers.
-}
