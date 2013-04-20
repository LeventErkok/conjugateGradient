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
-- applies to matrices that are symmetric and positive definite. See <http://en.wikipedia.org/wiki/Conjugate_gradient_method>
-- for details.
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
-- >>> import Data.IntMap.Strict
-- >>> import System.Random
-- >>> import Math.ConjugateGradient
-- >>> let a = SM (2, fromList [(0, SV (fromList [(0, 4), (1, 1)])), (1, SV (fromList [(0, 1), (1, 3)]))]) :: SM Double
-- >>> let b = SV (fromList [(0, 1), (1, 2)]) :: SV Double
-- >>> let g = mkStdGen 12345
-- >>> let x = solveCG g a b
-- >>> putStrLn $ showSolution 4 a b x
--       A       |   x    =   b   
-- --------------+----------------
-- 4.0000 1.0000 | 0.0909 = 1.0000
-- 1.0000 3.0000 | 0.6364 = 2.0000
---------------------------------------------------------------------------------

module Math.ConjugateGradient(
          -- * Types
            SV(..), SM(..)
          -- * Sparse operations
          , lookupSV, lookupSM, addSV, subSV, dotSV, normSV, sMulSV, sMulSM, mulSMV
          -- * Conjugate-Gradient solver
          , solveCG
          -- * Displaying solutions
          , showSolution
        ) where

import Data.List                          (intercalate)
import qualified Data.IntMap.Strict as IM (IntMap, lookup, map, unionWith, intersectionWith, fromList, findWithDefault, foldl')
import System.Random                      (Random, RandomGen, randomRs)
import Numeric                            (showFFloat)

-- | A sparse vector containing elements of type 'a'. Only the indices that contain non-@0@ elements should be given
-- for efficiency purposes. (Nothing will break if you put in elements that are @0@'s, it's just not as efficient.)
newtype SV a = SV (IM.IntMap a)

-- | A sparse matrix is essentially an int-map containing sparse row-vectors:
--
--     * The first element, @n@, is the number of rows in the matrix, including those with all @0@ elements.
--
--     * The matrix is implicitly assumed to be @nxn@, indexed by keys @(0, 0)@ to @(n-1, n-1)@.
--
--     * When constructing a sparse-matrix, only put in rows that have a non-@0@ element in them for efficiency.
--
--     * Note that you have to give all the non-0 elements: Even though the matrix must be symmetric for the algorithm
--       to work, the matrix should contain all the non-@0@ elements, not just the upper (or the lower)-triangle.
--
--     * Make sure the keys of the int-map is a subset of @[0 .. n-1]@, both for the row-indices and the indices of the vectors representing the sparse-rows.
newtype SM a = SM (Int, IM.IntMap (SV a))

---------------------------------------------------------------------------------
-- Sparse vector/matrix operations
---------------------------------------------------------------------------------

-- | Look-up a value in a sparse-vector.
lookupSV :: Num a => Int -> SV a -> a
lookupSV k (SV v) = IM.findWithDefault 0 k v

-- | Look-up a value in a sparse-matrix.
lookupSM :: Num a => (Int, Int) -> SM a -> a
lookupSM (i, j) (SM (_, m)) = maybe 0 (j `lookupSV`) (i `IM.lookup` m)

-- | Multiply a sparse-vector by a scalar.
sMulSV :: Num a => a -> SV a -> SV a
sMulSV s (SV v) = SV (IM.map (s *) v)

-- | Multiply a sparse-matrix by a scalar.
sMulSM :: Num a => a -> SM a -> SM a
sMulSM s (SM (n, m)) = SM (n, IM.map (s `sMulSV`) m)

-- | Add two sparse vectors.
addSV :: Num a => SV a -> SV a -> SV a
addSV (SV v1) (SV v2) = SV (IM.unionWith (+) v1 v2)

-- | Subtract two sparse vectors.
subSV :: Num a => SV a -> SV a -> SV a
subSV v1 (SV v2) = addSV v1 (SV (IM.map ((-1)*) v2))

-- | Dot product of two sparse vectors.
dotSV :: Num a => SV a -> SV a -> a
dotSV (SV v1) (SV v2) = IM.foldl' (+) 0 $ IM.intersectionWith (*) v1 v2

-- | Multiply a sparse matrix (nxn) with a sparse vector (nx1), obtaining a sparse vector (nx1).
mulSMV :: Num a => SM a -> SV a -> SV a
mulSMV (SM (_, m)) v = SV (IM.map (`dotSV` v) m)

-- | Norm of a sparse vector. (Square-root of its dot-product with itself.)
normSV :: RealFloat a => SV a -> a
normSV (SV v) = sqrt . IM.foldl' (\s e -> s + e*e) 0 $ v

-- | Conjugate Gradient Solver for the system @Ax=b@. See: <http://en.wikipedia.org/wiki/Conjugate_gradient_method>.
--
-- NB. Assumptions on the input:
--
--    * The @A@ matrix is symmetric and positive definite.
--
--    * All non-@0@ rows are present. (Even if the input is assumed symmetric, all rows must be present.)
--
--    * The indices start from @0@ and go consecutively up-to @n-1@. (Only non-@0@ value/row
--      indices has to be present, of course.)
--
-- For efficiency reasons, we do not check that these properties hold of the input. (If these assumptions are
-- violated, the algorithm will still produce a result, but not the one you expected!)
--
-- We perform either @10^6@ iterations of the Conjugate-Gradient algorithm, or until the error
-- factor is less than @1e-10@. The error factor is defined as the difference of the norm of
-- the current solution from the last one, as we go through the iterative solver. See
-- <http://en.wikipedia.org/wiki/Conjugate_gradient_method#Convergence_properties_of_the_conjugate_gradient_method>
-- for a discussion on the convergence properties of this algorithm.
--
-- The solver can throw an error if it does not converge by @10^6@ iterations. This is typically an indication
-- that the input matrix is not well formed, i.e., not symmetric positive-definite.
solveCG :: (RandomGen g, RealFloat a, Random a)
        => g          -- ^ The seed for the random-number generator.
        -> SM a       -- ^ The @A@ sparse matrix (@nxn@).
        -> SV a       -- ^ The @b@ sparse vector (@nx1@).
        -> SV a       -- ^ The @x@ sparse matrix (@nx1@), such that @Ax = b@.
solveCG g a@(SM (n, _)) b = cg a b x0
  where rs = take n (randomRs (0, 1) g)
        x0 = SV $ IM.fromList [p | p@(_, j) <- zip [0..] rs, j /= 0]

-- | The Conjugate-gradient algorithm. Our implementation closely follows the
-- one given here: <http://en.wikipedia.org/wiki/Conjugate_gradient_method#Example_code_in_Matlab>
cg :: RealFloat a => SM a -> SV a -> SV a -> SV a
cg a b x0 = cgIter (1000000 :: Int) (norm r0) r0 r0 x0
 where r0 = b `subSV` (a `mulSMV` x0)
       cgIter 0 _   _ _ _ = error "Conjugate Gradient: No convergence after 10^6 iterations. Make sure the input matrix is symmetric positive-definite!"
       cgIter i eps p r x
        -- Stop if the square of the error is less than 1e-20, i.e.,
        -- if the error itself is less than 1e-10.
        | eps' < 1e-20 = x'
        | True         = cgIter (i-1) eps' p' r' x'
        where ap    = a `mulSMV` p
              alpha = eps / (ap `dotSV` p)
              x'    = x `addSV` (alpha `sMulSV` p)
              r'    = r `subSV` (alpha `sMulSV` ap)
              eps'  = norm r'
              p'    = r' `addSV` ((eps' / eps) `sMulSV` p)
       norm (SV v) = IM.foldl' (\s e -> s + e*e) 0 v -- square of normSV, but no need for expensive square-root

-- | Display a solution in a human-readable form. Needless to say, only use this
-- method when the system is small enough to fit nicely on the screen.
showSolution :: RealFloat a
             => Int   -- ^ Precision: Use this many digits after the decimal point.
             -> SM a  -- ^ The @A@ matrix, @nxn@
             -> SV a  -- ^ The @b@ matrix, @nx1@
             -> SV a  -- ^ The @x@ matrix, @nx1@, as returned by 'solveCG', for instance.
             -> String
showSolution prec ma@(SM (n, _)) vb vx = intercalate "\n" $ header ++ res
  where res   = zipWith3 row a x b
        range = [0..n-1]
        sf d = showFFloat (Just prec) d ""
        a = [[sf ((i, j) `lookupSM` ma) | j <- range] | i <- range]
        x = [sf (i `lookupSV` vx) | i <- range]
        b = [sf (i `lookupSV` vb) | i <- range]
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

---------------------------------------------------------------------------------------------
-- Specialize for Float and Double instances
{-# SPECIALISE INLINE lookupSV :: Int        -> SV Float  -> Float                        #-}
{-# SPECIALISE INLINE lookupSV :: Int        -> SV Double -> Double                       #-}
{-# SPECIALISE INLINE lookupSM :: (Int, Int) -> SM Float  -> Float                        #-}
{-# SPECIALISE INLINE lookupSM :: (Int, Int) -> SM Double -> Double                       #-}
{-# SPECIALISE INLINE sMulSV   :: Float     -> SV Float  -> SV Float                      #-}
{-# SPECIALISE INLINE sMulSV   :: Double    -> SV Double -> SV Double                     #-}
{-# SPECIALISE INLINE sMulSM   :: Float     -> SM Float  -> SM Float                      #-}
{-# SPECIALISE INLINE sMulSM   :: Double    -> SM Double -> SM Double                     #-}
{-# SPECIALISE INLINE addSV    :: SV Float  -> SV Float  -> SV Float                      #-}
{-# SPECIALISE INLINE addSV    :: SV Double -> SV Double -> SV Double                     #-}
{-# SPECIALISE INLINE subSV    :: SV Float  -> SV Float  -> SV Float                      #-}
{-# SPECIALISE INLINE subSV    :: SV Double -> SV Double -> SV Double                     #-}
{-# SPECIALISE INLINE dotSV    :: SV Float  -> SV Float  -> Float                         #-}
{-# SPECIALISE INLINE dotSV    :: SV Double -> SV Double -> Double                        #-}
{-# SPECIALISE INLINE mulSMV   :: SM Float  -> SV Float  -> SV Float                      #-}
{-# SPECIALISE INLINE mulSMV   :: SM Double -> SV Double -> SV Double                     #-}
{-# SPECIALISE INLINE normSV   :: SV Float  -> Float                                      #-}
{-# SPECIALISE INLINE normSV   :: SV Double -> Double                                     #-}
{-# SPECIALISE        solveCG  :: RandomGen g => g -> SM Float  -> SV Float  -> SV Float  #-}
{-# SPECIALISE        solveCG  :: RandomGen g => g -> SM Double -> SV Double -> SV Double #-}
{-# SPECIALISE INLINE cg       :: SM Float  -> SV Float  -> SV Float  -> SV Float         #-}
{-# SPECIALISE INLINE cg       :: SM Double -> SV Double -> SV Double -> SV Double        #-}
