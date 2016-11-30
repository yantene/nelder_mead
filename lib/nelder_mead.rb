class Proc
  def nelder_mead(*init, accuracy)
    n = init.size
    alp = 1
    gam = 0.75 - 1.0 / (2 * n)
    rho = 1.0 + 2.0 / n
    sig = 1.0 - 1.0 / n

    nto = 20

    x = Array.new(n + 1){|i|
          init.map.with_index{|j, idx|
            j * (i - 1 == idx ? 1.2 : 1)
          }
        }

    f = Hash.new{|h, k| h[k] = self.call(*k)}

    lxo = []
    loop do
      x.sort_by!{|k| f[k]}

      xo = x[0..n-1].transpose.map{|xs| xs.inject(:+) / n}
      xr = xo.zip(x[n]).map{|xos, xns| xos + alp * (xos - xns)}

      if f[x[0]] <= f[xr] && f[xr] < f[x[n - 1]]
        # 反射
        x[n] = xr
      elsif f[xr] < f[x[0]]
        # 膨張
        xe = xo.zip(x[n]).map{|xos, xns| xos + gam * (xos - xns)}
        if f[xe] < f[xr]
          x[n] = xe
        else
          x[n] = xr
        end
      else
        # 収縮
        xc = xo.zip(x[n]).map{|xos, xns| xos + rho * (xos - xns)}
        if f[xc] < f[x[n]]
          x[n] = xc
        else
          1.upto(n) do |i|
            x[i] = x[0].zip(x[i]).map{|x0s, xis| x0s + sig * (xis - x0s)}
          end
        end
      end

      lxo = ([xo] + lxo).take(nto)
      break if lxo.size == nto && lxo.transpose.map{|k| k.each_cons(2).all?{|l| l.inject(:-).abs < accuracy}}.all?
    end
    x[0..n-1].transpose.map{|xs| xs.inject(:+) / n}
  end

  alias :downhill_simplex :nelder_mead
  alias :amoeba :nelder_mead
end
