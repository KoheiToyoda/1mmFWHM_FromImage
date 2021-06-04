using Base: Float64
using FileIO, Images, ImageView
using CSV, DataFrames
using Plots
using LsqFit
gr()

# 照射位置を決定する。
# #1 5x5 pixelの枠を作り、枠内最大強度となる点を掃引 
# 掃引範囲は画面右側とする(出射端の画面左端の方が強くなるので)
# #2 5x5の中心位置を入射面として、横方向1㎜(ピクセル換算する)して縦方向の強度分布をだす。
# #3 ガウスフィッティングしてFWHM出す。

#1 
function check_Intensity(img_path, peaktocheck)
    imgcolor = load(img_path)
    #imggray = Gray.(imgcolor)
    #greenやredでRGBパラメータ情報を返す
    img = green.(imgcolor)
    square = 5 
    size_ver, size_hor= size(img)
    start_hor = Int(size_hor/2)
    start_ver = 1
    end_hor = size_hor - square
    end_ver = size_ver - square
    @show start_hor,start_ver
    @show end_hor,end_ver
    ver_range = range(start_ver,stop = end_ver)
    hor_range = range(start_hor,stop = end_hor)
    
    max_sum_temp = 0
    max_x = 0
    max_y = 0

    for y in ver_range
        for x in hor_range
            # juliaにおける@view は指定した場所を参照する
            # 新たな行列を生成する。C言語におけるポインタのようなもの。
            # 実際にコピーした変数を用意するわけではない。
            # viewで作った行列を変更すると参照元も変更される
            # これを使うことでメモリ節約、高速化に貢献できる。

            clipped_img = @view img[y:y+square, x:x+square]
            sum_clipped_img = sum(clipped_img)
            #@show sum_clipped_img
            if max_sum_temp < sum_clipped_img
                max_sum_temp = sum_clipped_img
                max_x = x
                max_y = y
            end
        end
    end
    @show max_x,max_y

    # yの情報は後のフィッティングの時に使う。
    sweep_x = max_x - peaktocheck
    # 2
    intensity_array = @view img[:,sweep_x]

    xrange = range(1,stop = size_ver)
    pl1 = plot(xrange, intensity_array)
    # 関数途中のplotはdisplayで明示しないと表示されない。
    # 仮にdisplayを外すと210 にpl1 201にpl2の複合グラフになる。
    display(plot(pl1))
    # 光
    left = Base.prompt("山の左端を整数で入力")
    right = Base.prompt("山の右端を整数で入力")
    @show left, right
    # @view を使ってこういうこともできる。
    # (使わなくてもできる。)
    left_Array = @view intensity_array[begin:parse(Int,left)]
    right_Array = @view intensity_array[parse(Int,right):end]
    fill!(left_Array,0)
    fill!(right_Array,0)
    pl2 = plot(intensity_array)

    # LSQfit ここから最小二乗法で（と言ってもパッケージで）ガウスフィットをする。
    # model関数はガウシアンフィットを用意する。
    @. model(x,p) = p[3]*1/(p[1]*√(2π))*exp(-(x-p[2])^2/(2*(p[1])^2))
    # 初期値
    # pl はeltypeである必要(浮動小数型ってこと)
    # p0 = [σ, 中心位置, 振幅]
    p0 = [0.10, float(max_y), 1.0]
    fit = curve_fit(model, xrange , intensity_array, p0)
    @show fit.param
    # plot!は直前のPlot に上書きする。
    pl2 = plot!(xrange, model(range(1,stop = size_ver), fit.param))
    plot(pl2)
    # FWHM
    # http://hooktail.sub.jp/mathInPhys/fwhmsigma/
    return 2fit.param[1]*√(2log(2))
end

peaktocheck = 76
pix2mm = 1/76
 #76pixel = 1mm
file_base = "./1mmFWHM_FromImage/images/f=0."
file_branch = ["00","01","02","03","04","05","06"]
# file_branch = ["00"]
file_exte = ".jpeg"
time = [0,10,20,30,40,50,60]
FWHM_pix = Float64[]
FWHM_mm = Float64[]
for (index, bra) in enumerate(file_branch)
    img_path = file_base * bra * file_exte
    FWHM_now = check_Intensity(img_path, peaktocheck)
    push!(FWHM_pix,FWHM_now)
    push!(FWHM_mm ,FWHM_now*pix2mm)
end
@show time,FWHM_mm,FWHM_pix

#dataframe は辞書なしの方が楽
data = Dict(
    :time => time,
    :FWHM_mm => FWHM_mm,
    :FWHM_pix =>FWHM_pix,
)
df = DataFrame(data)
@show df

CSV.write("gauss_f0.csv", df) 