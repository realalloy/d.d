A new micro clearfix hack
The clearfix hack is a popular way to contain floats without resorting to using presentational markup. This article presents an update to the clearfix method that further reduces the amount of CSS required.

Demo: Micro clearfix hack

Known support: Firefox 3.5+, Safari 4+, Chrome, Opera 9+, IE 6+

The “micro clearfix” method is suitable for modern browsers and builds upon Thierry Koblentz’s “clearfix reloaded”, which introduced the use of both the :before and :after pseudo-elements.

Here is the updated code (I’ve used a shorter class name too):

/**
 * For modern browsers
 * 1. The space content is one way to avoid an Opera bug when the
 *    contenteditable attribute is included anywhere else in the document.
 *    Otherwise it causes space to appear at the top and bottom of elements
 *    that are clearfixed.
 * 2. The use of `table` rather than `block` is only necessary if using
 *    `:before` to contain the top-margins of child elements.
 */
.cf:before,
.cf:after {
    content: " "; /* 1 */
    display: table; /* 2 */
}

.cf:after {
    clear: both;
}

/**
 * For IE 6/7 only
 * Include this rule to trigger hasLayout and contain floats.
 */
.cf {
    *zoom: 1;
}
This “micro clearfix” generates pseudo-elements and sets their display to table. This creates an anonymous table-cell and a new block formatting context that means the :before pseudo-element prevents top-margin collapse. The :after pseudo-element is used to clear the floats. As a result, there is no need to hide any generated content and the total amount of code needed is reduced.

Including the :before selector is not necessary to clear the floats, but it prevents top-margins from collapsing in modern browsers. This has two benefits:

It ensures visual consistency with other float containment techniques that create a new block formatting context, e.g., overflow:hidden
It ensures visual consistency with IE 6/7 when zoom:1 is applied.
N.B.: There are circumstances in which IE 6/7 will not contain the bottom margins of floats within a new block formatting context. Further details can be found here: Better float containment in IE using CSS expressions.

The use of content:" " (note the space in the content string) avoids an Opera bug that creates space around clearfixed elements if the contenteditable attribute is also present somewhere in the HTML. Thanks to Sergio Cerrutti for spotting this fix. An alternative fix is to use font:0/0 a.

Legacy Firefox
Firefox < 3.5 will benefit from using Thierry’s method with the addition of visibility:hidden to hide the inserted character. This is because legacy versions of Firefox need content:"." to avoid extra space appearing between the body and its first child element, in certain circumstances (e.g., jsfiddle.net/necolas/K538S/.)

Alternative float-containment methods that create a new block formatting context, such as applying overflow:hidden or display:inline-block to the container element, will also avoid this behaviour in legacy versions of Firefox.
