// ==============================================================================
// GSEAlens Canvas Drag & Drop Logic
// ==============================================================================

$(document).ready(function() {

  // 监听父容器的事件
  $(document).on('mousedown', '.canvas-card', function(e) {

    // 只有点击头部手柄才能拖动
    if (!$(e.target).hasClass('card-header') && !$(e.target).hasClass('fa-arrows-alt')) {
      return;
    }

    e.preventDefault();

    var $card = $(this);
    var $canvas = $('#canvas_container');

    // 计算相对偏移量
    var offset = $card.offset();
    var rel_x = e.pageX - offset.left;
    var rel_y = e.pageY - offset.top;

    // 提升层级
    $card.css('z-index', 100);

    // 鼠标移动
    $(document).on('mousemove.canvas_drag', function(e) {
      var new_left = e.pageX - rel_x;
      var new_top = e.pageY - rel_y;

      // 边界检测 (可选)
      // ...

      $card.offset({ top: new_top, left: new_left });
    });

    // 鼠标松开
    $(document).on('mouseup.canvas_drag', function(e) {
      $(document).off('mousemove.canvas_drag');
      $(document).off('mouseup.canvas_drag');
      $card.css('z-index', 1);
    });
  });
});
