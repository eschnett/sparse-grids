{-# LANGUAGE TypeApplications #-}

{-# OPTIONS_GHC -Wno-type-defaults #-}

import qualified Test.Tasty
import Test.Tasty.Hspec

import qualified Data.Map as M
import SparseGrid



mkMap :: [([Int], a)] -> M.Map IndexL a
mkMap xs = M.fromList (map idx xs)
  where idx (is, x) = (idx_from_list is, x)

ord :: Order -> [([Int], AppLevel M.Map a)] -> AppOrder M.Map a
ord n ls = AppOrder n (mkMap ls)

lev :: [([Int], a)] -> AppLevel M.Map a
lev vs = AppLevel (mkMap vs)

mkAList :: [([Int], a)] -> AList IndexL a
mkAList xs = AList (map idx xs)
  where idx (is, x) = (idx_from_list is, x)

orda :: Order -> [([Int], AppLevel AList a)] -> AppOrder AList a
orda n ls = AppOrder n (mkAList ls)

leva :: [([Int], a)] -> AppLevel AList a
leva vs = AppLevel (mkAList vs)



main :: IO ()
main = do
    test <- testSpec "sparse-grids" spec
    Test.Tasty.defaultMain test

spec :: Spec
spec = parallel $ do
  it "td1" $ do
    td1 @Double `shouldBe`
      [ ord 1 []
      , ord 2 [([1,1], lev [([1,1],6.0)])]
      , ord 3 [ ([1,2], lev [([1,1],4.5), ([1,3],4.5)])
              , ([2,1], lev [([1,1],4.5), ([3,1],4.5)])
              ]
      ]
  it "take 3 $ make_discrete 1 fn0" $ do
    (take 3 $ make_discrete @Double 1 fn0) `shouldBe`
      [ ord 1 [([1], lev [([1],2.0)])]
      , ord 2 [([2], lev [([1],1.0), ([3],1.0)])]
      , ord 3 [([3], lev [([1],0.5), ([3],1.5), ([5],1.5), ([7],0.5)])]
      ]
  it "take 4 $ make_discrete 2 fn0" $ do
    (take 4 $ make_discrete @Double 2 fn0) `shouldBe`
      [ ord 1 []
      , ord 2 [([1,1], lev [([1,1],6.0)])]
      , ord 3 [ ([1,2], lev [([1,1],3.0), ([1,3],3.0)])
              , ([2,1], lev [([1,1],3.0), ([3,1],3.0)])
              ]
      , ord 4 [ ([1,3], lev [([1,1],1.5), ([1,3],4.5), ([1,5],4.5), ([1,7],1.5)])
              , ([2,2], lev [([1,1],1.5), ([1,3],1.5), ([3,1],1.5), ([3,3],1.5)])
              , ([3,1], lev [([1,1],1.5), ([3,1],4.5), ([5,1],4.5), ([7,1],1.5)])]
      ]
  it "eval_discrete_fn" $ do
    eval_discrete_fn (make_discrete @Double 2 fn0) (idx_from_list [2,2]) (idx_from_list [1,1]) `shouldBe` 1.5
  it "eval_discrete_fn" $ do
    eval_discrete_fn (make_discrete @Double 2 fn0) (idx_from_list [1,3]) (idx_from_list [1,5]) `shouldBe` 4.5
  it "eval_discrete_fn" $ do
    eval_discrete_fn (make_discrete @Float 2 fn0) (idx_from_list [1,3]) (idx_from_list [1,5]) `shouldBe` 4.5
  it "take 3 $ make_interpolant 1 fn0" $ do
    (take 3 $ make_interpolant @Double 1 fn0) `shouldBe`
      [ orda 1 [([1], leva [([1],2.0)])]
      , orda 2 [([2], leva [([1],0.0), ([3],0.0)])]
      , orda 3 [([3], leva [([1],0.0), ([3],0.0),([5],0.0), ([7],0.0)])]
      ]

-- TODO: Add all the tests from the main source code
