/**
 *  @file libcomp/geometry/circle.h
 */
#pragma once
#include "libcomp/geometry/point.hpp"

namespace lc {

/**
 *  @defgroup circle Circle
 *  @ingroup  geometry_primitives
 *  @{
 */

/**
 *  @brief 円
 */
struct Circle {
	/// 中心の座標
	Point c;
	/// 円の半径
	double r;

	/**
	 *  @brief コンストラクタ
	 *  @param[in] c  中心の座標
	 *  @param[in] r  円の半径
	 */
	explicit Circle(const Point &c = Point(), double r = 0.0) :
		c(c), r(r)
	{ }

	/**
	 *  @brief 不正な値を示す円の取得
	 *  @return 不正な円
	 */
	static Circle invalid(){
		double qnan = std::numeric_limits<double>::quiet_NaN();
		return Circle(Point::invalid(), qnan);
	}
	/**
	 *  @brief 円データが不正なものではないかの確認
	 *  @retval true   *thisが有効な円データであった場合
	 *  @retval false  *thisが無向な円データであった場合
	 */
	bool is_valid() const { return c.is_valid() && !std::isnan(r); }

	/**
	 *  @brief 円の比較 (<, 厳密評価)
	 *
	 *  コンテナで使用するためのもので数学的な意味はないことに注意。
	 *
	 *  @param[in] i      比較する値
	 *  @retval    true   *thisがiより辞書順で小さい場合
	 *  @retval    false  *thisがiより辞書順で大きい場合
	 */
	bool operator<(const Circle &i) const {
		return (r == i.r) ? (c < i.c) : (r < i.r);
	}

	/**
	 *  @brief 円と点の内外判定
	 *  @param[in] p   判定する点
	 *  @retval    0   p が *this の外側にある場合
	 *  @retval    1   p が *this の内側にある場合
	 *  @retval    -1  p が *this の円周上にある場合
	 */
	int contains(const Point &p) const {
		const double d2 = (p - c).norm(), r2 = r * r;
		if(d2 > r2){ return 0; }
		if(d2 < r2){ return 1; }
		return -1;
	}
};

/**
 *  @brief 円の比較 (==, 誤差許容)
 *  @param[in] a      比較する値
 *  @param[in] b      比較する値
 *  @retval    true   aとbが同じ円を表している場合
 *  @retval    false  aとbが同じ円を表している場合
 */
inline bool tolerant_eq(const Circle &a, const Circle &b){
	return tolerant_eq(a.c, b.c) && tolerant_eq(a.r, b.r);
}

/**
 *  @}
 */

}

