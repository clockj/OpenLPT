# OpenLPT

# Error Message
- 1: size error
- 2: type error
- 3: out of range
- 4: no enough space

# New version
- For Tsai's projection model: Xf,Yf start from 1
  But in C++, index starts from 0  
  Therefore, in my new code, I modify the projection model to make Xf,Yf start from 0   
