def mean(a)
  sum = 0
  for e in a do
    sum += e
  end
  return sum/a.length
end

def var(a)
  m = mean(a)
  var = 0
  for e in a do
    diff = m - e
    var += diff*diff
  end
  return (var/a.length)
end

# p "平均値"

# a = Array([0.1,0.2,0.3])

# p mean(a)

# a = Array([13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15])

# p mean(a)

a = Array([])
for i in 0..10000
  a.push(i)
end
p mean(a)

# p "分散"

# a = Array([0.1,0.2,0.3])

# p var(a)

# a = Array([13,14,7,12,10,6,8,15,4,14,9,6,10,12,5,12,8,8,12,15])

# p var(a)

# a = Array([])
# for i in 0..10000
#   a.push(i)
# end
# p var(a)
