// Add logo
remark.macros.addlogo = function (percentage) {
  var url = "img/logo.png";
  return '<img src="' + url + '" style="width: ' + percentage + '" />';
};


// Scale image (available on Remarkjs GH page)
remark.macros.scale = function (percentage) {
  var url = this;
  return '<img src="' + url + '" style="width: ' + percentage + '" />';
};


// toupper (see https://github.com/gnab/remark/issues/72)
remark.macros.upper = function () {
  // `this` is the value in the parenthesis, or undefined if left out
  return this.toUpperCase();
};


remark.macros.custom_hr = function(width, height) {
  return '<html><div style="float:left"></div><hr color="#e0b94c" style="margin-top:-60px" size=1px width=720px></html>'
}

