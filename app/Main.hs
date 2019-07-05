{-# LANGUAGE TypeApplications #-}

import qualified SparseGrid

main :: IO ()
main = do putStrLn "Sparse Grid interpolation"
          let (exact, appx, appxs) = SparseGrid.ta11 @Double
          putStrLn $ "exact solution: " ++ show exact
          putStrLn $ "best approximation: " ++ show appx
          putStrLn $ "successive contributions: " ++ show appxs
