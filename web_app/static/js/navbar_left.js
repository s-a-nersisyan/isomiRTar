$(window).on("load", function() {
  $('body').scrollspy({ target: '#navbar-left', offset: $("h2").first().offset().top + $(window).height()/2 });

  $(window).on('activate.bs.scrollspy', function () {
    history.replaceState({}, "", $('.nav-item .active').attr("href"));
  });
});
