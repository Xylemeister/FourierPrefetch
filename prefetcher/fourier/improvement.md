## Idea  1 (SMS)
1. we have a virtual page 4kb 
2. 64 cachelines 
3. for each page (or group of page) 
4. tehre is an origin 
5. we can use FFT to find/locate performant offsets (from the origin or line)
6. get all offsets and prefetch all of them


## Idea 2 (iterative)
1. strap fourier to SMS
2. simplify the patterns that we keep track of
3. questions on updating SMS entries 
4. IF we could use FFT to prune this update
5. DSPatch
