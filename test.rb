p "平均値"
p [0.1,0.2,0.3].mean

p [13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15].mean

a = Array([])
for i in 0..10000
  a.push(i)
end
p a.mean
p "分散"
p [0.1,0.2,0.3].var
p [13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15].var

a = Array([])
for i in 0..10000
  a.push(i)
end
p a.var
p "標準偏差"
p [1,2,3,4].std
p [0.1,0.2,0.3,0.4].std
p [13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15].std

p "kai mesure"
p [10,4].kai_mesure([2,4], 0.05)
p [160,40].kai_mesure([120,80], 0.05)
p [160,40].kai_mesure([120,80], 0.01)
p "t mesure"

a =[13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15] 

p a.t_mesure(12, 0.05)
p a.t_mesure(10, 0.05)

p "error case"
p ["aaa",2,3].mean
