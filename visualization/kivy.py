"""
CardioTRaNs 2019
Author: Rohit Suratekar

Kivy GUI functions
"""

from SecretColors import palette
from kivy.app import App
from kivy.lang import Builder
from kivy.properties import BooleanProperty
from kivy.uix.behaviors import FocusBehavior
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.recycleboxlayout import RecycleBoxLayout
from kivy.uix.recycleview import RecycleView
from kivy.uix.recycleview.layout import LayoutSelectionBehavior
from kivy.uix.recycleview.views import RecycleDataViewBehavior

p = palette.Palette(color_mode="rgb")

Builder.load_file('visualization/design.ky')


class SelectableRecycleBoxLayout(FocusBehavior,
                                 LayoutSelectionBehavior,
                                 RecycleBoxLayout):
    """
    Adds selection and focus behaviour to the view.
    """


class SelectableLabel(RecycleDataViewBehavior, Label):
    """
    Add selection support to the Label
    """
    index = None
    selected = BooleanProperty(False)
    selectable = BooleanProperty(True)

    r, g, b = p.gray(grade=10)
    row_color_normal = (r, g, b, 1)

    r, g, b = p.black()
    row_color_selected = (r, g, b, 1)

    def refresh_view_attrs(self, rv, index, data):
        """
        Catch and handle the view changes
        """
        self.index = index
        return super(SelectableLabel, self).refresh_view_attrs(
            rv, index, data)

    def on_touch_down(self, touch):
        """
        Add selection on touch down
        """
        if super(SelectableLabel, self).on_touch_down(touch):
            return True
        if self.collide_point(*touch.pos) and self.selectable:
            return self.parent.select_with_touch(self.index, touch)

    def apply_selection(self, rv, index, is_selected):
        """
        Respond to the selection of items in the view.
        """
        self.selected = is_selected
        if is_selected:
            rv.add_data(index)


class RV(RecycleView):
    selected_item = -1
    r, g, b = p.gray()
    base_color = (r, g, b, 1)

    def __init__(self, **kwargs):
        super(RV, self).__init__(**kwargs)
        for i in range(10):
            self.data.append({
                "text": "initial {}".format(i),
                "b_color": self.base_color}
            )

    def add_data(self, item: int):
        if self.selected_item != item:
            # self.data.append({"text": "now", "b_color": self.base_color})
            self.selected_item = item


class Test(BoxLayout):
    r, g, b = p.cool_gray(grade=10)
    box_background = (r, g, b, 1)

    def add_data(self):
        pass


class TestApp(App):

    def build(self):
        return Test()


def test():
    TestApp().run()


def run():
    test()
