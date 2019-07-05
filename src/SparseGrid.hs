{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE UndecidableInstances #-}

{-# OPTIONS_GHC -Wno-incomplete-uni-patterns #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}
{-# OPTIONS_GHC -Wno-type-defaults #-}

-- | This code was originally obtained on 2019-07-04 from
-- <http://okmij.org/ftp/Haskell/AlgorithmsH1.html>, which points to
-- <http://okmij.org/ftp/Haskell/SparseGrid.hs>. According to
-- <http://okmij.org/ftp/README.html>, the original is in the public
-- domain. I thank Oleg Kiselyov <oleg-at-okmij.org>.
--
-- For a background on sparse grids, see
--
-- * Hans-Joachim Bungartz and Michael Griebel: Sparse grids, Acta
--   Numerica (2004), pp. 1–123. DOI: 10.1017/S0962492904000182.
--   <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.9392&rep=rep1&type=pdf>.
--
-- * Jochen Garcke: Sparse Grid Tutorial. August 2006
--
--
--
-- Sparse hierarchical grid decomposition
-- for a function RI^d -\> R vanishing at the boundary.
-- Here RI is a closed interval on R.

module SparseGrid where

import Control.Monad
import Data.Bits (shiftL)
import qualified Data.Map as M
import qualified Data.Vector.Storable as S
import Foreign.Storable

default (Int)

-- | Aliases to make signatures more informative
type Dim = Int                          -- ^ d \in [0..Dim-1]
type Order = Int                        -- ^ Approximation order, n >=0

-- | Multi-dimensional index
newtype IDX = IDX{unIDX :: S.Vector Int}
  deriving (Eq, Ord)

instance Show IDX where
    show (IDX x) = show $ S.toList x

-- | Grid levels. There are 2^IndexL[i] grid points along the
-- dimension i. Each component of IndexL is >=1 [!]
type IndexL = IDX

-- | IndexJ[i] is the index of a grid point along the dimension i
-- (or the index of a basis function)
type IndexJ = IDX

-- | List of components for an IDX
idx_to_list :: IDX -> [Int]
idx_to_list = S.toList . unIDX

idx_from_list :: [Int] -> IDX
idx_from_list = IDX . S.fromList

-- | Compute the L1 norm of IDX
norm1 :: IDX -> Order
norm1 = S.sum . unIDX


choose :: [a] -> [a]
choose = id                             -- msum

-- | Compute all IndexL for the given Order.
-- That is, compute all l::IndexL such that norm1 l == n
levels_of_order :: Dim -> Order -> [IndexL]
levels_of_order d n = map idx_from_list $ makel d n
  where
  makel d n | n <= 0 || d <= 0 = mzero
  makel 1 n = return [n]
  makel d n = do
     li <- choose [1..n-1]
     ls <- makel (d-1) (n-li) 
     return $ li:ls

{-
tlevel1 = levels_of_order 2 3
*SparseGrid> levels_of_order 2 3
[[1,2],[2,1]]
*SparseGrid> levels_of_order 2 4
[[1,3],[2,2],[3,1]]
-}


-- | Within the code, we deal exclusively with functions defined on
-- the hyper-cube [0,1]^d. To remind us of that and avoid stupid
-- errors, we introduce the newtype DOMU. It has no run-time
-- overhead

type DOM a = S.Vector a            -- ^ unnormalized domain RI^d
newtype DOMU a = DOMU{domu :: S.Vector a} -- ^ [0,1]^d

{-
-- Convert the function RI^d -> R into [0,1]^d -> R

toDOMU :: [(a,a)] -> (DOM a -> a) -> (DOMU a -> a)
toDOMU minmax f = 

  x * span_i + min 
-}

-- | Approximant is a list of successive approximation orders
type Approximant coll a = [AppOrder coll a]

-- | One approximation of order n; a collection of levels IndexL
data AppOrder coll a = AppOrder !Order (coll IndexL (AppLevel coll a))

-- | One particular anisotropic grid Omega_L
-- A collection of points/basis functions identified by J and the values
-- at these points
newtype AppLevel coll a = AppLevel (coll IndexJ a)

deriving instance Eq (coll IndexL (AppLevel coll a)) => Eq (AppOrder coll a)
deriving instance Ord (coll IndexL (AppLevel coll a)) => Ord (AppOrder coll a)

deriving instance Eq (coll IndexJ a) => Eq (AppLevel coll a)
deriving instance Ord (coll IndexJ a) => Ord (AppLevel coll a)

instance Show a => Show (AppOrder M.Map a) where
    show (AppOrder _ coll) = "L" ++ show (M.toList coll)

instance Show a => Show (AppLevel M.Map a) where
    show (AppLevel coll) = "J" ++ show (M.toList coll)

newtype AList k v = AList{alist:: [(k,v)]}
  deriving (Eq, Ord)

instance Show a => Show (AppOrder AList a) where
    show (AppOrder _ coll) = "L" ++ show (alist coll)
instance Show a => Show (AppLevel AList a) where
    show (AppLevel coll) = "J" ++ show (alist coll)


-- ------------------------------------------------------------------------
-- * Computing the Interpolant
--
-- Representation of the function (DOMU a -> a) on the hierarchical,
-- discrete domain
-- A point in the discrete domain is identified by an index
-- (L,J) representing a fractional number J_i/2^L_i, i=0..Dim-1

make_discrete :: (Fractional a, Storable a)
              => Dim -> (DOMU a -> a) -> Approximant M.Map a
make_discrete d fn = map mk_order [1..]
 where
 mk_order n = AppOrder n (M.fromList . map mk_level $ levels_of_order d n)
 mk_level l = (l, AppLevel (M.fromList . map (mk_pt l) $ proper_js l))
 mk_pt l j = (j, fn $ pt_to_domu l j)
 pt_to_domu :: (Fractional a, Storable a) => IndexL -> IndexJ -> DOMU a
 pt_to_domu l j = DOMU $ S.zipWith dpoint_to_R (unIDX l) (unIDX j)


-- | Evaluate a discretized function at the given discrete point
eval_discrete_fn :: Approximant M.Map a -> IndexL -> IndexJ -> a
eval_discrete_fn grid l j = 
  let AppOrder _ ls = grid !! (norm1 l - 1)
      Just (AppLevel js) = M.lookup l ls
      Just r = M.lookup j js
  in r


-- | Given a level L, compute all proper points at that level
-- That is, compute J_i where J_i is odd and J_i < 2^L_i
proper_js :: IndexL -> [IndexJ]
proper_js ls = map idx_from_list . mapM mkj $ idx_to_list ls
 where
 mkj l = choose [1,3..pow2 l]

{-
*SparseGrid> proper_js (idx_from_list [1,2])
[[1,1],[1,3]]
*SparseGrid> proper_js (idx_from_list [2,3])
[[1,1],[1,3],[1,5],[1,7],[3,1],[3,3],[3,5],[3,7]]
-}

-- | Normalize a discrete point
-- Given (L,J), compute (L',J') that denotes the same Omega point
-- but each J'_i is odd
-- Return Nothing if the point is a boundary point (that is,
-- for some i, J_i = 0 or J_i = 2^L_i)
normalize_pt :: IndexL -> IndexJ -> Maybe (IndexL,IndexJ)
normalize_pt l j = do
   ljs <- zipWithM normalize (idx_to_list l) (idx_to_list j)
   let (ls,js) = unzip ljs
   return (idx_from_list ls, idx_from_list js)
 where
 normalize 0 _ = mzero
 normalize l j | (j',0) <- j `divMod` 2 = normalize (l-1) j'
 normalize l j = return (l,j)

-- | Given (not a boundary) point (L,J), compute the kernel
-- [-1/2 1 -1/2]^d centered at (L,J)
-- Return the array of normalized points and their coefficients
-- If the point falls at the boundary, ignore it.
neighbors_pt :: Fractional a => IndexL -> IndexJ -> [(a, (IndexL,IndexJ))]
neighbors_pt l j = do
 (f,js) <- mkk (idx_to_list j)
 maybe mzero (return . (,) f) $ normalize_pt l (idx_from_list js)
 where mkk [] = return (1,[])
       mkk (j:js) = do
          (f,j) <- choose [(-0.5,j-1), (1,j), (-0.5,j+1)]
          (fs,js) <- mkk js
          return (f*fs, j:js)

               
{-
*SparseGrid> take 3 $ make_interpolant 1 fn0
[L[([1],J[([1],2.0)])],L[([2],J[([1],0.0),([3],0.0)])],L[([3],J[([1],0.0),([3],0.0),([5],0.0),([7],0.0)])]]
*SparseGrid> neighbors_pt (idx_from_list [1,1]) (idx_from_list [1,1])
[(1.0,([1,1],[1,1]))]
*SparseGrid> neighbors_pt (idx_from_list [1,2]) (idx_from_list [1,1])
[(1.0,([1,2],[1,1])),(-0.5,([1,1],[1,1]))]
*SparseGrid> neighbors_pt (idx_from_list [2,2]) (idx_from_list [1,1])
[(1.0,([2,2],[1,1])),(-0.5,([2,1],[1,1])),(-0.5,([1,2],[1,1])),(0.25,([1,1],[1,1]))]
*SparseGrid> neighbors_pt (idx_from_list [2,2]) (idx_from_list [3,1])
[(-0.5,([1,2],[1,1])),(0.25,([1,1],[1,1])),(1.0,([2,2],[3,1])),(-0.5,([2,1],[3,1]))]
*SparseGrid> neighbors_pt (idx_from_list [3,2]) (idx_from_list [1,1])
[(1.0,([3,2],[1,1])),(-0.5,([3,1],[1,1])),(-0.5,([2,2],[1,1])),(0.25,([2,1],[1,1]))]
*SparseGrid> neighbors_pt (idx_from_list [3,3]) (idx_from_list [3,3])
[(0.25,([2,2],[1,1])),(-0.5,([2,3],[1,3])),(0.25,([2,1],[1,1])),(-0.5,([3,2],[3,1])),(1.0,([3,3],[3,3])),(-0.5,([3,1],[3,1])),(0.25,([1,2],[1,1])),(-0.5,([1,3],[1,3])),(0.25,([1,1],[1,1]))]

-}


-- | Compute the interpolant
make_interpolant :: (Fractional a, Storable a)
                 => Dim -> (DOMU a -> a) -> Approximant AList a
make_interpolant d fn = map mk_order [1..]
 where
 discrete_fn = make_discrete d fn
 eval_fn = eval_discrete_fn discrete_fn
 mk_order n = AppOrder n (AList . map mk_level $ levels_of_order d n)
 mk_level l = (l, AppLevel (AList . map (mk_pt l) $ proper_js l))
 mk_pt l j = (j, sum . map (\ (f,pt) -> f * uncurry eval_fn pt) $
                 neighbors_pt l j)

-- ------------------------------------------------------------------------
-- | Evaluating the Approximant
--
-- Evaluate a basis function at level L
-- Given a point x in DOMU, compute the index J of the basis function
-- phi and the value of the J-th function at point x
-- In the dimension i:
--  the index j = 1 + 2*floor(2^(l-1) * x)
--  the function is 1 - abs(2*y - 1), y = fraction(2^(l-1) * x)

eval_phi :: (RealFrac a, Storable a) => IndexL -> DOMU a -> (IndexJ, a)
eval_phi l (DOMU x) = (idx_from_list js,r)
 where
 (js,r) = foldr (\ (j,v) (js,vs) -> (j:js,v * vs)) ([],1.0) $ 
          zipWith eval1 (idx_to_list l) (S.toList x)
 eval1 l x = let xscaled = scale_x (l-1) x
                 (xw,y)  = properFraction xscaled
                 j       = 1 + 2 * xw
             in (j,1 - abs(2*y-1))

{-
*SparseGrid> eval_phi (idx_from_list [1]) (to_DOMU [0.4])
([1],0.8)
*SparseGrid> eval_phi (idx_from_list [1]) (to_DOMU [0.5])
([1],1.0)
*SparseGrid> eval_phi (idx_from_list [2]) (to_DOMU [0.5])
([3],0.0)
*SparseGrid> eval_phi (idx_from_list [2]) (to_DOMU [0.6])
([3],0.3999999999999999)
*SparseGrid> eval_phi (idx_from_list [2]) (to_DOMU [0.4])
([1],0.3999999999999999)
*SparseGrid> eval_phi (idx_from_list [2]) (to_DOMU [0.25])
([1],1.0)
*SparseGrid> eval_phi (idx_from_list [3]) (to_DOMU [0.25])
([3],0.0)
*SparseGrid> eval_phi (idx_from_list [1,1]) (to_DOMU [0.4,0.4])
([1,1],0.6400000000000001)
*SparseGrid> eval_phi (idx_from_list [1,2]) (to_DOMU [0.4,0.4])
([1,1],0.31999999999999995)
*SparseGrid> eval_phi (idx_from_list [1,3]) (to_DOMU [0.4,0.4])
([1,3],0.6399999999999999)
*SparseGrid> eval_phi (idx_from_list [1,3]) (to_DOMU [0.4,0.25])
([1,3],0.0)
*SparseGrid> eval_phi (idx_from_list [1,3]) (to_DOMU [0.5,0.4])
([1,3],0.7999999999999998)
-}

-- | Evaluate an approximant of the given order at the given point in DOMU
eval_approximant_order :: (RealFrac a, Storable a)
                       => AppOrder AList a -> DOMU a -> a
eval_approximant_order (AppOrder _ (AList ls)) xs = foldr eval_l 0 ls
 where
 eval_l (l,AppLevel (AList js)) acc = 
     let (j,phiv) = eval_phi l xs
         (Just a) = lookup j js
     in acc + a * phiv

-- | Here, we assume the grid list to be finite
-- We should truncate the approximant first
eval_approximant :: (RealFrac a, Storable a)
                 => Approximant AList a -> DOMU a -> a
eval_approximant grid xs = 
    foldr (\l acc -> acc + eval_approximant_order l xs) 0 grid


-- ------------------------------------------------------------------------
-- * Utilities

-- | Compute 2^l. We assume that Order is <= 31 
pow2 :: Order -> Int
pow2 l = 1 `shiftL` l

-- | compute j/2^l
dpoint_to_R :: Fractional a => Int -> Int -> a
dpoint_to_R l j = fromIntegral j / 2^l

-- | Scale a point: x * 2^l
scale_x :: Fractional a => Int -> a -> a
scale_x l x = x * 2^l

-- | Used for testing. Assume [a] is in [0,1]^d
to_DOMU :: Storable a => [a] -> DOMU a
to_DOMU = DOMU . S.fromList

-- ------------------------------------------------------------------------
-- * Tests

-- | Sample function on a 2d DOMU

fn0 :: (Num a, Storable a) => DOMU a -> a
fn0 (DOMU dm) = product $ 
                zipWith (\c x -> c * (1 - abs (2*x - 1))) [2,3] (S.toList dm)

fn1 :: (Num a, Storable a) => DOMU a -> a
fn1 (DOMU dm) = product $ 
                zipWith (\c x -> c * (1 - (2*x - 1)^2)) [2,3] (S.toList dm)

td1 :: (Fractional a, Storable a) => [AppOrder M.Map a]
td1 = take 3 $ make_discrete 2 fn1
{-
 [L[], -- order 1
  L[([1,1],J[([1,1],6.0)])], -- order 2
  -- order 3
  L[([1,2],J[([1,1],4.5),([1,3],4.5)]),
    ([2,1],J[([1,1],4.5),([3,1],4.5)])]]
-}

{-
*SparseGrid> take 3 $ make_discrete 1 fn0
 [L[([1],J[([1],2.0)])],
  L[([2],J[([1],1.0),([3],1.0)])],
  L[([3],J[([1],0.5),([3],1.5),([5],1.5),([7],0.5)])]]


*SparseGrid> take 4 $ make_discrete 2 fn0
 [L[],
  L[([1,1],J[([1,1],6.0)])],
  L[([1,2],J[([1,1],3.0),([1,3],3.0)]),([2,1],J[([1,1],3.0),([3,1],3.0)])],
  L[([1,3],J[([1,1],1.5),([1,3],4.5),([1,5],4.5),([1,7],1.5)]),
    ([2,2],J[([1,1],1.5),([1,3],1.5),([3,1],1.5),([3,3],1.5)]),
    ([3,1],J[([1,1],1.5),([3,1],4.5),([5,1],4.5),([7,1],1.5)])]]

-- at (1.4,1/4)
*SparseGrid> eval_discrete_fn (make_discrete 2 fn0) (idx_from_list [2,2]) (idx_from_list [1,1])
1.5
*SparseGrid> eval_discrete_fn (make_discrete 2 fn0) (idx_from_list [1,3]) (idx_from_list [1,5])
4.5

-}

{-
*SparseGrid> take 3 $ make_interpolant 1 fn0
[L[([1],J[([1],2.0)])], -- level 1
 L[([2],J[([1],0.0),([3],0.0)])], -- level 2, etc.
 L[([3],J[([1],0.0),([3],0.0),([5],0.0),([7],0.0)])]]
-}

ta0 :: (Fractional a, Storable a) => [AppOrder AList a]
ta0 = take 4 $ make_interpolant 2 fn0
{-
 [L[], -- order 1
  L[([1,1],J[([1,1],6.0)])], -- order 2
  -- order 3
  L[([1,2],J[([1,1],0.0),([1,3],0.0)]),([2,1],J[([1,1],0.0),([3,1],0.0)])],
  L[([1,3],J[([1,1],0.0),([1,3],0.0),([1,5],0.0),([1,7],0.0)]),
    ([2,2],J[([1,1],0.0),([1,3],0.0),([3,1],0.0),([3,3],0.0)]),
    ([3,1],J[([1,1],0.0),([3,1],0.0),([5,1],0.0),([7,1],0.0)])]]
-}

ta01 :: (RealFrac a, Storable a) => a
ta01 = eval_approximant ta0 (to_DOMU [0.5,0.5])
-- 6.0

ta02 :: (RealFrac a, Storable a) => a
ta02 = eval_approximant ta0 (to_DOMU [0.49,0.51])
-- 5.7623999999999995

ta1 :: (Fractional a, Storable a) => [AppOrder AList a]
ta1 = take 4 $ make_interpolant 2 fn1
{-
 [L[], -- order 1
  L[([1,1],J[([1,1],6.0)])], -- order 2
  -- order 3
  L[([1,2],J[([1,1],1.5),([1,3],1.5)]),([2,1],J[([1,1],1.5),([3,1],1.5)])],
  -- order 4
  L[([1,3],J[([1,1],0.375),([1,3],0.375),([1,5],0.375),([1,7],0.375)]),
    ([2,2],J[([1,1],0.375),([1,3],0.375),([3,1],0.375),([3,3],0.375)]),
    ([3,1],J[([1,1],0.375),([3,1],0.375),([5,1],0.375),([7,1],0.375)])]]
-}

-- | Show successive contributions
eval_app :: (RealFrac a, Storable a)
         => [AppOrder AList a] -> DOMU a -> [a]
eval_app grid xs = map (\l -> eval_approximant_order l xs) grid

ta11 :: RealFrac a => Storable a => (a, a, [a])
ta11 = let xs = to_DOMU [0.49,0.51]
           exact = fn1 xs
           -- successive approximations
           appxs = take 9 $ eval_app (make_interpolant 2 fn1) xs
           appx = sum appxs
       in (exact, appx, appxs)
{-
(5.995200960000001,
 5.994952734375,
 [0.0,5.7623999999999995,0.1176000000000001,5.940000000000005e-2,3.0000000000000023e-2,1.515000000000001e-2,7.650000000000006e-3,2.254687499999997e-3,4.980468750000017e-4])
-}
