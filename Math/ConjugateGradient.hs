-- http://en.wikipedia.org/wiki/Conjugate_gradient_method
module Math.ConjugateGradient(M, V, solve, dispSolution) where

import Data.Maybe (fromMaybe)
import qualified Data.IntMap as IM
import System.Random
import Data.List (intercalate)
import Numeric

type FP = Double

dispSolution :: Int -> Int -> M -> V -> V -> String
dispSolution padLen n ma vb vx = intercalate "\n" $ zipWith3 row a b x
  where a = [[ma `mLookup` (i, j) | j <- [0..n-1]] | i <- [0..n-1]]
        b = [vb `vLookup` i | i <- [0..n-1]]
        x = [vx `vLookup` i | i <- [0..n-1]]
        row :: [FP] -> FP -> FP -> String
        row as bv xv = unwords (map sh as) ++ " | " ++ sh bv ++ " = " ++ sh xv
        sh :: FP -> String
        sh d = pad $ showFFloat (Just 3) d ""
        pad s = reverse $ take (length s `max` padLen) $ reverse s ++ repeat ' '

type V = IM.IntMap FP
type M = IM.IntMap V

mLookup :: M -> (Int, Int) -> FP
mLookup m (i, j) = maybe 0 (`vLookup` j) (i `IM.lookup` m)

vLookup :: V -> Int -> FP
vLookup m k = fromMaybe 0 (k `IM.lookup` m)

-- a: nxn
-- b: nx1
-- x: nx1
solve :: RandomGen g => g -> Int -> M -> V -> V
solve g n a b = cg a b x0
  where rs = take n (randomRs (0, 1) g)
        x0 = IM.fromList [p | p@(_, j) <- zip [0..] rs, j /= 0]

mulMV :: M -> V -> V
mulMV m v = IM.fromList [(i, val) | i <- IM.keys m, let val = mul i, val /= 0]
  where mul i  = maybe (error "mulMV: Impossible happened!") elt (i `IM.lookup` m)
        elt cv = sum [d * v `vLookup` j | (j, d) <- IM.toList cv]

subV :: V -> V -> V
subV v1 v2 = addV v1 (IM.map ((-1)*) v2)

addV :: V -> V -> V
addV = IM.unionWith (+)

dot :: V -> V -> FP
dot v1 v2 = IM.fold (+) 0 $ IM.intersectionWith (*) v1 v2

norm :: V -> FP
norm = IM.fold (\e s -> e*e + s) 0

sMulV :: FP -> V -> V
sMulV s = IM.map (s *)

cg :: M -> V -> V -> V
cg a b x0 = cgIter 1001 (norm r0) r0 r0 x0
 where r0 = subV b (mulMV a x0)
       cgIter :: Int -> FP -> V -> V -> V -> V
       cgIter 0 _   _ _ x = x
       cgIter i eps p r x
        | sqrt eps' < 1e-8 = x'
        | True             = cgIter (i-1) eps' p' r' x'
        where ap    = a `mulMV` p
              alpha = eps / dot ap p
              x'    = x `addV` (alpha `sMulV` p)
              r'    = r `subV` (alpha `sMulV` ap)
              eps'  = norm r'
              p'    = r' `addV` ((eps' / eps) `sMulV` p)
